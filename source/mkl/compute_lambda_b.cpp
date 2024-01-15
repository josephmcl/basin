
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
    da.mode = SPARSE_FILL_MODE_UPPER;
	da.diag = SPARSE_DIAG_NON_UNIT;

    std::size_t findex, gindex;
    std::vector<std::tuple<std::size_t, std::size_t>> ind;
    for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
        for (std::size_t j = 0; j != sbp.n_blocks; ++j) {
            findex = FT_symbols[i][j];
            gindex = (i * sbp.n) + j;
            if (findex != 0) {
                ind.push_back({findex - 1, gindex});
            }
        }
    }
    
    real_t *mg, *lb;
    std::size_t j, k;
    for (std::size_t i = 0; i != ind.size(); ++i) {
        j = std::get<0>(ind[i]);
        k = std::get<1>(ind[i]);
        lb = &λb[k];
        mg = &Mg[k];
        status = mkl_sparse_d_mv(
            SPARSE_OPERATION_NON_TRANSPOSE, -1., Fsparse[j], da,
            mg, 
            1., lb);
        mkl_sparse_status(status);
    }
}
