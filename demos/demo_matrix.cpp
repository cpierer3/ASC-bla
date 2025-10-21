#include <iostream>
#include <matrix.hpp>
#include <cassert>

namespace bla = ASC_bla;

void test_row_major_and_col_major() {
  bla::Matrix<double, bla::RowMajor> A(3, 3);
  bla::Matrix<double, bla::ColMajor> B(3, 2);
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      A(i, j) = i * 3 + j;
    }
  }
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 2; j++) {
      B(i, j) = i * 2 + j;
    }
  }

  std::cout << "Row Major Matrix A:\n" << A << std::endl;
  std::cout << "Column Major Matrix B:\n" << B << std::endl;

  std::cout << "A: ";
  for (size_t i =0; i<3*3; i++) {
    std::cout << A.getRawDataDanger()[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << "B: ";
  for (size_t i =0; i<3*2; i++) {
    std::cout << B.getRawDataDanger()[i] << ", ";
  }
  std::cout << std::endl;

  // get rows
  for (size_t i = 0; i < 3; i++) {
    auto rowA = A.row(i);
    auto rowB = B.row(i);
    std::cout << "Row " << i << " of A: " << rowA << std::endl;
    std::cout << "Row " << i << " of B: " << rowB << std::endl;
  }
  // get cols
  for (size_t j = 0; j < 3; j++) {
    auto colA = A.col(j);
    std::cout << "Column " << j << " of A: " << colA << std::endl;
    if (j >= 2) continue;
    auto colB = B.col(j);
    std::cout << "Column " << j << " of B: " << colB << std::endl;
  }
}

void test_multiplications() {
  using namespace std;
  namespace bla = ASC_bla;

  // A is 2x3
  bla::Matrix<double, bla::RowMajor> A(2,3);
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;

  // B is 3x2
  bla::Matrix<double, bla::RowMajor> B(3,2);
  B(0,0) = 7;  B(0,1) = 8;
  B(1,0) = 9;  B(1,1) = 10;
  B(2,0) = 11; B(2,1) = 12;

  cout << "A:\n" << A << endl;
  cout << "B:\n" << B << endl;
  cout << "a rows: " << A.rows() << ", cols: " << A.cols() << endl;
  cout << "b rows: " << B.rows() << ", cols: " << B.cols() << endl;

  bla::Matrix<double, bla::RowMajor> C = A * B; // expected 2x2: [[58,64],[139,154]]
  cout << "A*B =\n" << C << endl;
  assert(C.rows() == 2 && C.cols() == 2);
  assert(std::abs(C(0,0) - 58.0) < 1e-12);
  assert(std::abs(C(0,1) - 64.0) < 1e-12);
  assert(std::abs(C(1,0) - 139.0) < 1e-12);
  assert(std::abs(C(1,1) - 154.0) < 1e-12);

  // vector multiplication: v has size 3
  bla::Vector<double> v(3);
  v(0) = 1; v(1) = 0; v(2) = -1;

  cout << "v = " << v << endl;
  auto Av = A * v; // expected size 2: [-2, -2]
  cout << "A*v = " << Av << endl;
  assert(Av.size() == 2);
  assert(std::abs(Av(0) + 2.0) < 1e-12);
  assert(std::abs(Av(1) + 2.0) < 1e-12);

  cout << "Matrix-matrix and matrix-vector multiplication tests passed.\n";

  auto A_scaled = 2.0 * A;
  assert(A_scaled(0,0) == 2.0);
  assert(A_scaled(1,2) == 12.0);

  cout << "Matrix-scalar tests passed.\n";

}

int main() {
  test_row_major_and_col_major();
  test_multiplications();
  size_t n = 3;
  size_t m = 3;
  bla::Matrix<double, bla::RowMajor> A(n, m);
  bla::Matrix<double, bla::RowMajor> B(n, m);
  bla::Matrix<double, bla::ColMajor> BC(n, m);

  for (size_t i = 0; i < A.rows(); i++) {
    for (size_t j = 0; j < A.cols(); j++) {
      A(i, j) = j + 1;
      B(i, j) = (i + 1) * (j + 1) * (j + 1);
      BC(i, j) = (i + 1) * (j + 1) * (j + 1);
    }
  }

  auto C = A + B;
  bla::Matrix<double, bla::RowMajor> D = C;
  std::cout << "D = A + B = " << D << std::endl;
}
