#include "poisson_2d.h"

#include "csr.h"
#include <iostream>

int main(int argc, char **argv) {

    // batch gemm. 
    //oneapi::mkl::blas::row_major::gemm_batch
    // cblas_dgemm_batch;

  /* Initialize libraries needed for simulation. */

    poisson_2d::problem(4,3);

    /*
    csr<int> x;

    x = csr<int>{4,6};
    
    x(7, 2, 4);
    x(6, 2, 3);
    x(8, 3, 5);

    x(2, 0, 1);
    x(3, 1, 1);
    x(4, 1, 3);
    x(5, 2, 2);
    x(1, 0, 0);
    

    for (auto &e: x.v) {
      std::cout << e << " ";
    } std::cout << std::endl;

    for (auto &e: x.c) {
      std::cout << e << " ";
    } std::cout << std::endl;
    for (auto &e: x.r) {
      std::cout << e << " ";
    } std::cout << std::endl;

    std::cout << x.nnz() << std::endl;
        */



}