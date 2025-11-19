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

#include <matrixview.hpp>
#include <mat_simd.hpp>
#include <lapack_interface.hpp>


int main() {

  std::srand(static_cast<unsigned int>(std::time(nullptr)));


  std::vector<std::pair<std::string, std::vector<std::pair<int, double>>>> all_data;
  std::vector<std::pair<std::string, std::pair<ASC_bla::ORDERING, ASC_bla::ORDERING>>> configs = {
    {"RowMajor x RowMajor", {ASC_bla::RowMajor, ASC_bla::RowMajor}},
    {"RowMajor x ColMajor", {ASC_bla::RowMajor, ASC_bla::ColMajor}},
    {"ColMajor x RowMajor", {ASC_bla::ColMajor, ASC_bla::RowMajor}},
    {"ColMajor x ColMajor", {ASC_bla::ColMajor, ASC_bla::ColMajor}}
  };

  for (int n = 2; n <= 2048; n *= 2) {
    // Generate random data once for this n
    std::vector<double> dataA(n * n), dataB(n * n);
    for (int i = 0; i < n * n; ++i) {
      dataA[i] = static_cast<double>(rand()) / RAND_MAX;
      dataB[i] = static_cast<double>(rand()) / RAND_MAX;
    }
    for (const auto& config : configs) {
      std::string label = config.first;
      auto A_ord = config.second.first;
      auto B_ord = config.second.second;
      std::vector<std::pair<int, double>>* data_ptr = nullptr;
      // Find or create the data vector for this config
      auto it = std::find_if(all_data.begin(), all_data.end(), [&](const auto& p) { return p.first == label; });
      if (it == all_data.end()) {
        all_data.push_back({label, {}});
        data_ptr = &all_data.back().second;
      } else {
        data_ptr = &it->second;
      }
      // Create matrices with the same data but different ordering
      if (A_ord == ASC_bla::RowMajor && B_ord == ASC_bla::RowMajor) {
        ASC_bla::Matrix<double, ASC_bla::RowMajor> A(n, n), B(n, n);
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j) {
            A(i, j) = dataA[i * n + j];
            B(i, j) = dataB[i * n + j];
          }
        double n_cubed = static_cast<double>(n) * n * n;
        int runs = 1 + static_cast<int>(std::min(1e8 / n_cubed, 1000.0));
        std::cout << "\n=== " << label << " n=" << n << " runs: " << runs << std::endl;
        ASC_bla::Matrix<double, ASC_bla::RowMajor> C(n, n);
        auto ts = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < runs; ++i) {
          C = Multi(A, B);
        }
        auto te = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = te - ts;
        double time_per_run = elapsed_seconds.count() / runs;
        std::cout << "n = " << std::setw(4) << n
                  << ", runs = " << std::setw(5) << runs
                  << ", time = " << std::fixed << std::setprecision(10) << time_per_run << " s"
                  << std::endl;
        data_ptr->push_back({n, time_per_run});
      } else if (A_ord == ASC_bla::RowMajor && B_ord == ASC_bla::ColMajor) {
        ASC_bla::Matrix<double, ASC_bla::RowMajor> A(n, n);
        ASC_bla::Matrix<double, ASC_bla::ColMajor> B(n, n);
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j) {
            A(i, j) = dataA[i * n + j];
            B(i, j) = dataB[i * n + j];
          }
        double n_cubed = static_cast<double>(n) * n * n;
        int runs = 1 + static_cast<int>(std::min(1e8 / n_cubed, 1000.0));
        std::cout << "\n=== " << label << " n=" << n << " runs: " << runs << std::endl;
        ASC_bla::Matrix<double, ASC_bla::RowMajor> C(n, n);
        auto ts = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < runs; ++i) {
          C = Multi(A, B);
        }
        auto te = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = te - ts;
        double time_per_run = elapsed_seconds.count() / runs;
        std::cout << "n = " << std::setw(4) << n
                  << ", runs = " << std::setw(5) << runs
                  << ", time = " << std::fixed << std::setprecision(10) << time_per_run << " s"
                  << std::endl;
        data_ptr->push_back({n, time_per_run});
      } else if (A_ord == ASC_bla::ColMajor && B_ord == ASC_bla::RowMajor) {
        ASC_bla::Matrix<double, ASC_bla::ColMajor> A(n, n);
        ASC_bla::Matrix<double, ASC_bla::RowMajor> B(n, n);
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j) {
            A(i, j) = dataA[i * n + j];
            B(i, j) = dataB[i * n + j];
          }
        double n_cubed = static_cast<double>(n) * n * n;
        int runs = 1 + static_cast<int>(std::min(1e8 / n_cubed, 1000.0));
        std::cout << "\n=== " << label << " n=" << n << " runs: " << runs << std::endl;
        ASC_bla::Matrix<double, ASC_bla::RowMajor> C(n, n);
        auto ts = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < runs; ++i) {
          C = Multi(A, B);
        }
        auto te = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = te - ts;
        double time_per_run = elapsed_seconds.count() / runs;
        std::cout << "n = " << std::setw(4) << n
                  << ", runs = " << std::setw(5) << runs
                  << ", time = " << std::fixed << std::setprecision(10) << time_per_run << " s"
                  << std::endl;
        data_ptr->push_back({n, time_per_run});
      } else if (A_ord == ASC_bla::ColMajor && B_ord == ASC_bla::ColMajor) {
        ASC_bla::Matrix<double, ASC_bla::ColMajor> A(n, n), B(n, n);
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j) {
            A(i, j) = dataA[i * n + j];
            B(i, j) = dataB[i * n + j];
          }
        double n_cubed = static_cast<double>(n) * n * n;
        int runs = 1 + static_cast<int>(std::min(1e8 / n_cubed, 1000.0));
        std::cout << "\n=== " << label << " n=" << n << " runs: " << runs << std::endl;
        ASC_bla::Matrix<double, ASC_bla::RowMajor> C(n, n);
        auto ts = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < runs; ++i) {
          C = Multi(A, B);
        }
        auto te = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = te - ts;
        double time_per_run = elapsed_seconds.count() / runs;
        std::cout << "n = " << std::setw(4) << n
                  << ", runs = " << std::setw(5) << runs
                  << ", time = " << std::fixed << std::setprecision(10) << time_per_run << " s"
                  << std::endl;
        data_ptr->push_back({n, time_per_run});
      }
    }
  }

  // Output all results
  for (const auto& entry : all_data) {
    std::cout << "\nResults for " << entry.first << ":\n";
    for (size_t i = 0; i < entry.second.size(); ++i) {
      std::cout << entry.second[i].first << ", " << std::fixed << std::setprecision(10) << entry.second[i].second << std::endl;
      if (i < entry.second.size() - 1) {
        std::cout << ", ";
      }
    }
  }

  return 0;
}
