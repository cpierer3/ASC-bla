# Welcome to ASC-bla's documentation!


ASC-bla is a C++ library for basic linear algebra operations.
The library provides template classes **Vector** and **Matrix**.

## Installation

install it via git-clone:

    git clone https://github.com/TUWien-ASC/ASC-bla.git


To configure and build some tests do

    cd ASC-bla
    mkdir build
    cd build
    cmake ..
    make
    

## Using ASC-bla

To use ASC-bla in your code, set the compiler include path properly, and include the header files

    #include <vector.hpp>
    #include <matrix.hpp>

All objects are implemented in the namespace ASC_bla. To use them with less typing, you can set

    namespace bla = ASC_bla;

or even

    
    using namespace ASC_bla;

    

You can create vectors and compute with vectors like:

                 
```cpp
Vector<double> x(5), y(5), z(5);
for (int i = 0; i < x.Size(); i++)
   x(i) = i;
y = 5.0
z = x+3*y;
cout << "z = " << z << endl;
```

For matrices you can choose between row-major (`RowMajor`) or column-major (`ColMajor`) storage,
default is row-major.

```cpp
Matrix<double,RowMajor> m1(5,3), m2(3,3);
for (int i = 0; i < m1.Height(); i++)
  for (int j = 0; j < m1.Width(); j++)
    m1(i,j) = i+j;
m2 = 3.7;
Matrix product = m1 * m2;
```

You can extract a row or a column from a matrix:

```cpp
Vector col1 = product.Col(1);
```

---

## Advanced Features

### Expression Templates for Efficient Operations

ASC-bla uses expression templates to enable efficient and elegant matrix and vector operations. This means you can write code like:

```cpp
Matrix<double> A(3,3), B(3,3), C(3,3);
// ... initialize A and B ...
C = A + B;         // Matrix addition
C = A - B;         // Matrix subtraction
C = 2.0 * A;       // Scalar multiplication
C = A * B;         // Matrix-matrix multiplication
```

These operations are performed efficiently without unnecessary temporary objects, thanks to the expression template mechanism in the library.

### Matrix Views: Submatrices, Rows, Columns, and Transpose

ASC-bla provides powerful matrix views, allowing you to work with submatrices, rows, columns, and transposed matrices without copying data:

```cpp
Matrix<double> M(5,5);
// ... initialize M ...

// Extract a row or column as a view
auto row2 = M.row(2);      // 3rd row as a vector view
auto col3 = M.col(3);      // 4th column as a vector view

// Extract a submatrix (rows 1 to 3)
auto subM = M.rows(1,4);   // rows 1,2,3 (end index exclusive)

// Extract a range of columns
auto subCols = M.cols(2,5); // columns 2,3,4

// Transpose view (no copy)
auto Mt = M.transpose();

// Assign to a view
row2 = 0.0;                // set all elements in row 2 to zero
subM = 1.0;                // set all elements in submatrix to one
```

### Fundamental Matrix operations such as inversion

You can compute the inverse of a square matrix using:

```cpp
Matrix<double> A(3,3);
// ... initialize A ...
Matrix<double> Ainv = A.inverse();
```

This uses Gauss-Jordan elimination with partial pivoting.

```

---

For more details, see the source files `matrixexpr.hpp` and `matrixview.hpp`.


