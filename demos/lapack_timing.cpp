//
// Created by jonas on 10/22/25.
//

#include <iostream>
#include <chrono>   // for high-res timing
#include <matrixview.hpp>
#include <lapack_interface.hpp>

void lapack(long n) {
  size_t flops = n * n * n;
  size_t runs = size_t(1e9 / flops) + 1;
  std::cout << "runs: " << runs << std::endl;

  auto start = std::chrono::high_resolution_clock::now();
  ASC_bla::Matrix<double, ASC_bla::ColMajor> A(n, n), B(n, n), C(n, n);
  for (size_t i = 0; i < runs; i++) {
    ASC_bla::multMatMatLapack(A, B, C);
  }
  auto end = std::chrono::high_resolution_clock::now();
  double time = std::chrono::duration<double>(end - start).count();

  std::cout << "lapack n = " << n << ", time = " << time << " s, GFlops = " << (flops * runs) / time * 1e-9 << std::endl;
}
void our(long n) {
  size_t flops = n * n * n;
  size_t runs = size_t(1e9 / flops) + 1;

  std::cout << "runs: " << runs << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  ASC_bla::Matrix<double, ASC_bla::ColMajor> A(n, n), B(n, n), C(n, n);
  for (size_t i = 0; i < runs; i++) {
    C = A * B;
  }
  auto end = std::chrono::high_resolution_clock::now();
  double time = std::chrono::duration<double>(end - start).count();

  std::cout << "our n = " << n << ", time = " << time << " s, GFlops = " << (flops * runs) / time * 1e-9 << std::endl;
}

int main() {
  int n = 1000;
  our(n);
  lapack(n);
}


