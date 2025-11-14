# Require C++17 for constexpr
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(LAPACK REQUIRED)

add_executable (demo_vector demo_vector.cpp)
target_sources (demo_vector PUBLIC ../src/vector.hpp ../src/vecexpr.hpp)

add_executable (demo_matrix demo_matrix.cpp)
target_sources (demo_matrix PUBLIC ../src/vector.hpp ../src/vecexpr.hpp ../src/matrixview.hpp ../src/matrixexpr.hpp)
add_executable (timing_demo timing_demo.cpp)
target_sources (timing_demo PUBLIC ../src/vector.hpp ../src/vecexpr.hpp ../src/matrixview.hpp ../src/matrixexpr.hpp)
target_link_libraries (timing_demo PUBLIC LAPACK::LAPACK)

target_compile_options(demo_matrix PRIVATE -save-temps -g0)

add_executable (matrix_timing matrix_timing.cpp)
target_compile_options(matrix_timing PRIVATE -save-temps -g0)
