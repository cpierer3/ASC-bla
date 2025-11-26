#ifndef FILE_MATSIMD
#define FILE_MATSIMD

#include <matrixview.hpp>
#include "../myASC-HPC/src/simd.hpp"
#include "../myASC-HPC/src/taskmanager.hpp"

namespace ASC_bla {
  template<typename T>
  Matrix<T, ORDERING::RowMajor>
  Multi(const MatrixView<T, ORDERING::RowMajor> A, const MatrixView<T, ORDERING::RowMajor> B) {
    Matrix<T, ORDERING::RowMajor> C(A.rows(), B.cols());

    constexpr int COLS = 16;
    constexpr int ROWS = 8;

    ASC_HPC::RunParallel(A.rows() / ROWS, [&C, &B, &A](int idx, int size) {
      const int i = idx * ROWS;
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
    });
    return C;
  }

}

#endif