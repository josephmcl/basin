
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

    std::size_t findex, j_global;
    double * lb, *mg;
    #pragma omp parallel for collapse(2) private(findex, lb, mg, status) num_threads(7)
    for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
        for (std::size_t j = 0; j != sbp.rank_limit_u; ++j) {
            j_global = sbp.rank_index_u[j];
            if (FT_symbols[i][j_global] > 0) {
                findex = FT_symbols[i][j_global] - 1;
                lb = &λb[i * sbp.n];
                // NOTE: j is the local (rank) index. j_global is the 
                //       global element index. we use j below because 
                //       Mg is a rank computation. But we j_global 
                //       elsewhere to look up where the global index 
                //       in other the symbolic data structures.
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
