# Require C++17 for constexpr
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(LAPACK REQUIRED)

add_executable(demo_vector demo_vector.cpp)
target_sources(demo_vector PUBLIC ../src/vector.hpp ../src/vecexpr.hpp)

add_executable(demo_matrix demo_matrix.cpp)
target_sources(demo_matrix PUBLIC ../src/vector.hpp ../src/vecexpr.hpp ../src/matrixview.hpp ../src/matrixexpr.hpp)
add_executable(lapack_timing lapack_timing.cpp)
target_sources(lapack_timing PUBLIC ../src/vector.hpp ../src/vecexpr.hpp ../src/matrixview.hpp ../src/matrixexpr.hpp)
target_link_libraries(lapack_timing PUBLIC LAPACK::LAPACK)

target_compile_options(demo_matrix PRIVATE -save-temps -g0)

add_executable(matrix_timing matrix_timing.cpp)
target_sources(matrix_timing PUBLIC
    ../myASC-HPC/concurrentqueue/concurrentqueue.h
    ../myASC-HPC/src/taskmanager.cpp
    ../myASC-HPC/src/timer.cpp ../src/vector.hpp
)
target_compile_options(matrix_timing PRIVATE -save-temps)
target_link_libraries(matrix_timing PUBLIC LAPACK::LAPACK)
