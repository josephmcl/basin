#include "compute_u.h"

void compute_u(
    real_t *u,
    vv<sparse_matrix_t> &M,
    real_t *rhs, 
    components &sbp) {  

    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

    sparse_status_t status;
    double *up, *rp;
    std::size_t k;
    std::size_t limit = sbp.n_blocks;   
    
    #pragma omp parallel for private(up, rp, k) num_threads(sbp.n_threads)
    for (std::size_t i = 0; i != limit; ++i) {
        auto td = omp_get_thread_num();
        up = &u[i * sbp.n * sbp.n];
        rp = &rhs[i * sbp.n * sbp.n];
        k = mi[i % sbp.n_blocks_dim];
        status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, up , sbp.n * sbp.n, 
            rp, sbp.n * sbp.n);
        mkl_sparse_status(status);
    }    
}