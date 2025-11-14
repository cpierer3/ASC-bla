#ifndef FILE_MATSIMD
#define FILE_MATSIMD

#include <matrixview.hpp>
#include "../myASC-HPC/src/simd.hpp"

namespace ASC_bla {
    template <typename t>
    Matrix Multi(const MatrixView<T> A, const MatrixView<T> B) {

        for (size_t i =0; i < A.rows()/2; i+= 2 ){
            A.row(i)
            ;


        }
        
    }

}

#endif