//
// Created by jonas on 11/12/25.
//
#include <iostream>
#include <vector>
#include <chrono>      // For high-precision timing
#include <cmath>       // For std::pow and std::min
#include <algorithm>   // For std::min
#include <iomanip>     // For std::setw, std::fixed, std::setprecision
#include <utility>     // For std::pair
#include <cstdlib>     // For rand()
#include <ctime>       // For time() to seed rand()

// Include your existing Matrix class header
// We'll use the "Matrix.h" stub provided below for a complete, runnable example.
#include <matrixview.hpp>
#include <mat_simd.hpp>
#include <lapack_interface.hpp>


int main() {
  // Seed the random number generator.
  // We initialize matrices with random data to prevent the
  // compiler from optimizing away the multiplication.
  std::srand(static_cast<unsigned int>(std::time(nullptr)));

  // This will store the (n, time) pairs, just like the Python 'data' list
  std::vector<std::pair<int, double>> data;

  std::cout << "Starting matrix multiplication benchmark..." << std::endl;
  std::cout << "--------------------------------------------------------" << std::endl;

  // The Python loop logic was:
  // n = 1 -> (loop) -> n = 2 (benchmark)
  // n = 2 -> (loop) -> n = 4 (benchmark)
  // ...
  // n = 1024 -> (loop) -> n = 2048 (benchmark)
  // This 'for' loop is the direct C++ equivalent.
  for (int n = 16; n <= 2048; n *= 2) {
//  for (int n = 8; n <= 8; n *= 2) {

    ASC_bla::Matrix<double, ASC_bla::RowMajor> A(n, n);
    ASC_bla::Matrix<double, ASC_bla::RowMajor> B(n, n);

    // Initialize matrices with random values
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
        A(i, j) = static_cast<double>(i+j+2.2);
        B(i, j) = static_cast<double>(i+j+1.1);
      }
    }

    // Calculate the number of runs
    // We use double for the calculation to avoid integer overflow on n*n*n
    double n_cubed = static_cast<double>(n) * n * n;
    int runs = 1 + static_cast<int>(std::min(1e8 / n_cubed, 1000.0));

    std::cout << "runs: " << runs << std::endl;
    // Pre-allocate the result matrix so its allocation isn't part of the timed loop
    ASC_bla::Matrix<double, ASC_bla::RowMajor> C(n, n);

    // Start the timer
    auto ts = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < runs; ++i) {
      C = Multi(A, B);
    }

    // Stop the timer
    auto te = std::chrono::high_resolution_clock::now();

    ASC_bla::Matrix<double, ASC_bla::RowMajor> D(n, n);
//    D = A * B; // Using overloaded operator*
    ASC_bla::multMatMatLapack(A, B, D);
    // Verify that C and D are approximately equal
    double max_diff = 0.1;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double diff = std::abs(C(i, j) - D(i, j));
        if (diff > max_diff) {
          std::cerr << "Difference at (" << i << ", " << j << "): C = " << C(i, j) << ", D = " << D(i, j) << std::endl;
        }
      }
    }

    // Calculate the total elapsed time in seconds
    std::chrono::duration<double> elapsed_seconds = te - ts;

    // Calculate the average time per run
    double time_per_run = elapsed_seconds.count() / runs;

    // Print the result for this 'n'
    std::cout << "n = " << std::setw(4) << n
              << ", runs = " << std::setw(5) << runs
              << ", time = " << std::fixed << std::setprecision(10) << time_per_run << " s"
              << std::endl;

    // Store the result
    data.push_back({n, time_per_run});
  }

  // Print the final data vector, mimicking the Python list output
  for (size_t i = 0; i < data.size(); ++i) {
    std::cout << data[i].first << ", " << std::fixed << std::setprecision(10) << data[i].second << std::endl;
    if (i < data.size() - 1) {
      std::cout << ", ";
    }
  }

  return 0;
}
