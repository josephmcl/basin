
#include "compute_lambda_b.h"

void compute_lambda_b(
    real_t *λb, 
    std::vector<sparse_matrix_t> &Fsparse, 
    real_t *Mg, 
    vv<std::size_t> &FT_symbols, 
    components &sbp) {

    sparse_status_t status;

    matrix_descr da;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;

    std::size_t findex, gindex;
    double * lb, *mg;
    #pragma omp parallel for collapse(2) private(findex, gindex, lb, mg, status) num_threads(sbp.n_threads)
    for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
        for (std::size_t j = 0; j != sbp.n_blocks; ++j) {
            if (FT_symbols[i][j] > 0) {
                findex = FT_symbols[i][j] - 1;
                lb = &λb[i * sbp.n];
                mg = &Mg[j * sbp.n * sbp.n];
                status = mkl_sparse_d_mv(
                    SPARSE_OPERATION_NON_TRANSPOSE, -1., 
                    Fsparse[findex], da,
                    mg, 
                    1., lb);
                mkl_sparse_status(status);
            }
        }
    }
}
