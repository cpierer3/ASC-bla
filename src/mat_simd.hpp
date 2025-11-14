#ifndef FILE_MATSIMD
#define FILE_MATSIMD

#include <matrixview.hpp>
#include "../myASC-HPC/src/simd.hpp"

namespace ASC_bla {
  template<typename T>
  Matrix<T, ORDERING::RowMajor> Multi(const MatrixView<T, ORDERING::RowMajor> A, const MatrixView<T, ORDERING::RowMajor> B) {
    Matrix<T, ORDERING::RowMajor> C(A.rows(), B.cols());


    constexpr int STRIDE = 8;
    for (size_t i = 0; i < A.rows() / 2; i += 2) {
      auto r1 = A.row(i);
      auto r2 = A.row(i + 1);
      for (size_t j = 0; j < B.cols() / STRIDE; j += STRIDE) {

        ASC_HPC::SIMD<double, STRIDE> a0 = ASC_HPC::SIMD<double, STRIDE>(0.0);
        ASC_HPC::SIMD<double, STRIDE> a1 = ASC_HPC::SIMD<double, STRIDE>(0.0);

        for (size_t k=0; k<A.cols(); k++) {
          auto c = ASC_HPC::SIMD<double, STRIDE>(B.data()[k*B.dist()+j]);
          a0 += r1(k) * c;
          a1 += r2(k) * c;
        }

        a0.store(&C.data()[i*C.dist()+j]);
        a1.store(&C.data()[(i+1)*C.dist()+j]);
      }
    }
    return C;
  }

}

#endif