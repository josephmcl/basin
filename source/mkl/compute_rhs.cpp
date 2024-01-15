#include "compute_rhs.h"

void compute_rhs(
    real_t *rhs,
    std::vector<sparse_matrix_t> &F,
    real_t *lambda,
    vv<std::size_t> &F_symbols,
    components &sbp) {

    // Compute rhs = - F * lambda

    sparse_status_t status;
    matrix_descr da;
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> temp;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    for (std::size_t i = 0; i != sbp.n_blocks; ++i) {
        for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
            if (F_symbols[i][j] != 0) {
                temp.push_back({i, j, F_symbols[i][j] - 1});
            }
        }
    }

    /*
    sparse_index_base_t indexing;
    MKL_INT rows, cols, *rowst, *rowe, *coli, *ia;
    real_t *vals;
    status = mkl_sparse_d_export_csr(
      F[0], &indexing, &rows, &cols, &rowst, &rowe, &coli, &vals);
    mkl_sparse_status(status);

    //ia = (MKL_INT *) MKL_malloc(sizeof(MKL_INT) * ((16) + 1), 64);
    //std::memcpy(&ia[0], &rowst[0], sizeof(MKL_INT) * 16);
    //ia[16] = rowe[16 - 1];

    std::cout << rows << " " << cols << std::endl;
    */

    std::size_t i, j, k;
    double *l, *r;
    for (std::size_t a = 0; a != temp.size(); ++a) {
        i = std::get<0>(temp[a]);
        j = std::get<1>(temp[a]);
        k = std::get<2>(temp[a]);
        l = &lambda[j * sbp.n];
        r = &rhs[i * sbp.n * sbp.n];
        status = mkl_sparse_d_mv(
            SPARSE_OPERATION_TRANSPOSE, -1., 
            F[k], da,
            l, 
            1., r);
        mkl_sparse_status(status); 
    }
} 