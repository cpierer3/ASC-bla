#ifndef FILE_MATSIMD
#define FILE_MATSIMD

#include <matrixview.hpp>
#include "../myASC-HPC/src/simd.hpp"

namespace ASC_bla {
  template<typename T>
  Matrix<T, ORDERING::RowMajor> Multi(const MatrixView<T, ORDERING::RowMajor> A, const MatrixView<T, ORDERING::RowMajor> B) {
    Matrix<T, ORDERING::RowMajor> C(A.rows(), B.cols());


    constexpr int COLS = 8;
    constexpr int ROWS = 4;
    for (size_t i = 0; i < A.rows(); i += ROWS) {
      for (size_t j = 0; j < B.cols(); j += COLS) {

        ASC_HPC::SIMD<double, COLS> a[ROWS];
        for (int t=0; t<ROWS; t++) {
          a[t] = ASC_HPC::SIMD<double, COLS>(0.0);
        }
        for (size_t k=0; k<A.cols(); k++) {
          auto c = ASC_HPC::SIMD<double, COLS>(&B.data()[k * B.dist() + j]);
          for (int t=0; t<ROWS; t++) {
            a[t] += A(i + t, k) * c;
          }
        }

        for (int t=0; t<ROWS; t++) {
          a[t].store(&C.data()[(i + t) * C.dist() + j]);
        }
      }
    }
    return C;
  }

  template<typename T>
  // naive testing for time difference
  Matrix <T, ORDERING::RowMajor> Multi(const MatrixView<T, ORDERING::ColMajor> A, const MatrixView<T, ORDERING::RowMajor> B) {
    //std::cout << "Using transpose multiplication!" << std::endl;
    return Multi(A.transpose(), B);
}

template<typename T>
  Matrix <T, ORDERING::RowMajor> Multi(const MatrixView<T, ORDERING::RowMajor> A, const MatrixView<T, ORDERING::ColMajor> B) {
    //std::cout << "Using transpose multiplication!" << std::endl;
    return Multi(A, B.transpose());
}

template<typename T>
  Matrix <T, ORDERING::RowMajor> Multi(const MatrixView<T, ORDERING::ColMajor> A, const MatrixView<T, ORDERING::ColMajor> B) {
    //std::cout << "Using transpose multiplication!" << std::endl;
    return Multi(A.transpose(), B.transpose());

}
}
#endif