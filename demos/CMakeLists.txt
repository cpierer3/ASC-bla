# Require C++17 for constexpr
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable (demo_vector demo_vector.cpp)
target_sources (demo_vector PUBLIC ../src/vector.hpp ../src/vecexpr.hpp)

add_executable (demo_matrix demo_matrix.cpp)
target_sources (demo_matrix PUBLIC ../src/vector.hpp ../src/vecexpr.hpp ../src/matrix.hpp ../src/matrixexpr.hpp)