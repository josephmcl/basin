#include "infrastructure.h"
#include "poisson_2d.h"

int main(int argc, char **argv) {

    // batch gemm. 
    //oneapi::mkl::blas::row_major::gemm_batch
    // cblas_dgemm_batch;

  /* Initialize libraries needed for simulation. */
  if (!infrastructure::initialize(argc, argv)) exit(-1);; { 

    poisson_2d::problem();

  // Cleanup petsc 
  } infrastructure::cleanup(); return 0;

}