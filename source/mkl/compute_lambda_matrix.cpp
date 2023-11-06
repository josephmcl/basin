#include "compute_lambda_matrix.h"

auto poisson_2d::compute_lambda_matrix(
        double     &*lambdaA,
  const double     **F,
  const double     **MinvFT,
  const std::size_t n_batch,
  const components &sbp) -> void {

    std::vector m(n_batch), n(n_batch), k(n_batch);
    for (std::size_t i = 0; i != n_batch; ++i) {
        m[i] = sbp.n * sbp.n;
        n[i] = sbp.n;
        k[i] = sbp.n;
    }


    void cblas_dgemm_batch (
        const CBLAS_LAYOUT Layout, 
        const CBLAS_TRANSPOSE* transa_array, 
        const CBLAS_TRANSPOSE* transb_array, 
        const MKL_INT* m_array, 
        const MKL_INT* n_array, 
        const MKL_INT* k_array, 
        const double* alpha_array, 
        const double **a_array, 
        const MKL_INT* lda_array, 
        const double **b_array, 
        const MKL_INT* ldb_array, 
        const double* beta_array, 
        double **c_array, 
        const MKL_INT* ldc_array, 
        const MKL_INT group_count, 
        const MKL_INT* group_size);

    /* mkl batch gemm. */
    cblas_dgemm_batch(
        CblasNoTrans,
        t_array,
        t_array,
        



    return;
}