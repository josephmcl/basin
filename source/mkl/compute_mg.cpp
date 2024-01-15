#include "compute_mg.h"

void compute_mg(
    real_t *Mg,
    std::vector<sparse_matrix_t> &M,
    real_t *g,
    components &sbp) {  


    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

    sparse_status_t status;
    double *gp, *mgp;
    std::size_t k;
    std::size_t limit = sbp.n_blocks;   
    for (std::size_t i = 0; i != limit; ++i) {
        gp = &g[i * sbp.n * sbp.n];
        mgp = &Mg[i * sbp.n * sbp.n];
        k = mi[i % sbp.n_blocks_dim];
        status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
        mkl_sparse_status(status);
    }
}