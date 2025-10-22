#include <iostream>
#include <cassert>
#include <matrixview.hpp>
#include <cmath>

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

void test_rows_view() {
  namespace bla = ASC_bla;

  bla::Matrix<double, bla::RowMajor> A(4, 3);
  bla::Matrix<double, bla::ColMajor> B(4, 3);
  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;
  A(2,0) = 7; A(2,1) = 8; A(2,2) = 9;
  A(3,0) = 10; A(3,1) = 11; A(3,2) = 12;
  B(0,0) = 1; B(0,1) = 2; B(0,2) = 3;
  B(1,0) = 4; B(1,1) = 5; B(1,2) = 6;
  B(2,0) = 7; B(2,1) = 8; B(2,2) = 9;
  B(3,0) = 10; B(3,1) = 11; B(3,2) = 12;
  auto A_rows_1_3 = A.rows(1, 3); // should be rows 1 and 2
  assert(A_rows_1_3.rows() == 2 && A_rows_1_3.cols() == 3);
  assert(A_rows_1_3(0,0) == 4 && A_rows_1_3(1,2) == 9);
  auto B_rows_0_2 = B.rows(1, 3); // should be rows 0 and 1
  assert(B_rows_0_2.rows() == 2 && B_rows_0_2.cols() == 3);
  assert(B_rows_0_2(0,0) == 4 && B_rows_0_2(1,2) == 9);
}

void test_cols_view() {
  namespace bla = ASC_bla;

  // Create a 3x4 matrix for testing column views
  bla::Matrix<double, bla::RowMajor> A(3, 4);
  bla::Matrix<double, bla::ColMajor> B(3, 4);

  // Fill matrices with distinct values
  A(0,0) = 1;  A(0,1) = 2;  A(0,2) = 3;  A(0,3) = 4;
  A(1,0) = 5;  A(1,1) = 6;  A(1,2) = 7;  A(1,3) = 8;
  A(2,0) = 9;  A(2,1) = 10; A(2,2) = 11; A(2,3) = 12;

  B(0,0) = 1;  B(0,1) = 2;  B(0,2) = 3;  B(0,3) = 4;
  B(1,0) = 5;  B(1,1) = 6;  B(1,2) = 7;  B(1,3) = 8;
  B(2,0) = 9;  B(2,1) = 10; B(2,2) = 11; B(2,3) = 12;

  // Test cols(1,3) - should get columns 1 and 2
  auto A_cols_1_3 = A.cols(1, 3);
  assert(A_cols_1_3.rows() == 3 && A_cols_1_3.cols() == 2);
  assert(A_cols_1_3(0,0) == 2 && A_cols_1_3(0,1) == 3);
  assert(A_cols_1_3(1,0) == 6 && A_cols_1_3(1,1) == 7);
  assert(A_cols_1_3(2,0) == 10 && A_cols_1_3(2,1) == 11);

  auto B_cols_1_3 = B.cols(1, 3);
  assert(B_cols_1_3.rows() == 3 && B_cols_1_3.cols() == 2);
  assert(B_cols_1_3(0,0) == 2 && B_cols_1_3(0,1) == 3);
  assert(B_cols_1_3(1,0) == 6 && B_cols_1_3(1,1) == 7);
  assert(B_cols_1_3(2,0) == 10 && B_cols_1_3(2,1) == 11);

  std::cout << "cols(1,3) view test passed.\n";
}

void test_transpose() {
  using namespace std;
  namespace bla = ASC_bla;

  // Create a 3x4 RowMajor matrix
  bla::Matrix<double, bla::RowMajor> A(3, 4);
  A(0,0) = 1;  A(0,1) = 2;  A(0,2) = 3;  A(0,3) = 4;
  A(1,0) = 5;  A(1,1) = 6;  A(1,2) = 7;  A(1,3) = 8;
  A(2,0) = 9;  A(2,1) = 10; A(2,2) = 11; A(2,3) = 12;

  auto At = A.transpose();

//  cout << "A (3x4):\n" << A << endl;
//  cout << "A^T (4x3):\n" << At << endl;
//
  // Verify dimensions are swapped
  assert(At.rows() == 4 && At.cols() == 3);

  // Verify elements: At(i,j) should equal A(j,i)
  for (size_t i = 0; i < At.rows(); ++i)
    for (size_t j = 0; j < At.cols(); ++j)
      assert(At(i,j) == A(j,i));

  // Test with ColMajor matrix
  bla::Matrix<double, bla::ColMajor> B(2, 3);
  B(0,0) = 1; B(0,1) = 2; B(0,2) = 3;
  B(1,0) = 4; B(1,1) = 5; B(1,2) = 6;

  auto Bt = B.transpose();

//  cout << "B (2x3):\n" << B << endl;
//  cout << "B^T (3x2):\n" << Bt << endl;

  assert(Bt.rows() == 3 && Bt.cols() == 2);
  for (size_t i = 0; i < Bt.rows(); ++i)
    for (size_t j = 0; j < Bt.cols(); ++j)
      assert(Bt(i,j) == B(j,i));

  cout << "Transpose test passed.\n";
}

void test_inverse() {
  using namespace std;
  namespace bla = ASC_bla;

  // Test 2x2 matrix inversion
  bla::Matrix<double, bla::RowMajor> A(2, 2);
  A(0,0) = 4;  A(0,1) = 7;
  A(1,0) = 2;  A(1,1) = 6;

  ASC_bla::Matrix<double, ASC_bla::RowMajor> Ainv = A.inverse();

  cout << "A:\n" << A << endl;
  cout << "A^(-1):\n" << Ainv << endl;

  // Verify A * A^(-1) = I
  auto I = A * Ainv;
  cout << "A * A^(-1):\n" << I << endl;

  assert(abs(I(0,0) - 1.0) < 1e-10);
  assert(abs(I(0,1)) < 1e-10);
  assert(abs(I(1,0)) < 1e-10);
  assert(abs(I(1,1) - 1.0) < 1e-10);

  // Test 3x3 matrix inversion
  bla::Matrix<double, bla::RowMajor> B(3, 3);
  B(0,0) = 1;  B(0,1) = 2;  B(0,2) = 3;
  B(1,0) = 0;  B(1,1) = 1;  B(1,2) = 4;
  B(2,0) = 5;  B(2,1) = 6;  B(2,2) = 0;

  auto Binv = B.inverse();

  cout << "B:\n" << B << endl;
  cout << "B^(-1):\n" << Binv << endl;

  // Verify B * B^(-1) = I
  auto I2 = B * Binv;
  cout << "B * B^(-1):\n" << I2 << endl;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      double expected = (i == j) ? 1.0 : 0.0;
      assert(abs(I2(i,j) - expected) < 1e-9);
    }
  }

  // Test with ColMajor matrix
  bla::Matrix<double, bla::ColMajor> C(2, 2);
  C(0,0) = 3;  C(0,1) = 2;
  C(1,0) = 1;  C(1,1) = 4;

  auto Cinv = C.inverse();
  auto I3 = C * Cinv;

  assert(abs(I3(0,0) - 1.0) < 1e-10);
  assert(abs(I3(0,1)) < 1e-10);
  assert(abs(I3(1,0)) < 1e-10);
  assert(abs(I3(1,1) - 1.0) < 1e-10);

  cout << "Inverse test passed.\n";
}

int main() {
  test_row_major_and_col_major();
  test_rows_view();
  test_cols_view();
  test_transpose();
  test_inverse();
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
