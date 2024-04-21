#include "compute_mg.h"

void compute_mg(
    real_t *Mg,
    std::vector<double *> &M,
    std::vector<int *> &Mpiv,
    real_t *g,
    components &sbp) {  

    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

    sparse_status_t status;
    double *gp, *mgp;
    std::size_t k;
    std::size_t limit = sbp.rank_limit_u;   
    // #pragma omp parallel for private(gp, mgp, k) num_threads(sbp.n_threads)
    for (std::size_t i = 0; i != limit; ++i) {
        auto td = omp_get_thread_num();
        gp = &g[i * sbp.n * sbp.n];
        mgp = &Mg[i * sbp.n * sbp.n];
        std::memcpy(mgp, gp, sbp.n * sbp.n * sizeof(double));

        // rank_index_u is the global index, thats what we modulo. 
        k = mi[sbp.rank_index_u[i] % sbp.n_blocks_dim];
        
        /*
        status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
        mkl_sparse_status(status);
        */
       
        auto res = LAPACKE_dgetrs(
          LAPACK_COL_MAJOR, 
          'N', // No transpose
          sbp.n * sbp.n, // Rows in A & b
          1, // Number of RHS
          M[k], 
          sbp.n * sbp.n, // lda 
          Mpiv[k], 
          mgp, // RHS
          sbp.n * sbp.n); // 

        if (res != 0) {
          std::cout << res << std::endl;
        }
    }
}


void compute_mg_mpi(
    real_t *Mg,
    vv<sparse_matrix_t> &M,
    real_t *g,
    components &sbp) {

    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;


    sparse_status_t status;
    double *gp, *mgp;
    std::size_t k;
    std::size_t limit = sbp.n_blocks;   
    #pragma omp parallel for private(gp, mgp, k) num_threads(sbp.n_threads)
    for (std::size_t i = 0; i != limit; ++i) {
        auto td = omp_get_thread_num();
        gp = &g[i * sbp.n * sbp.n];
        mgp = &Mg[i * sbp.n * sbp.n];
        k = mi[i % sbp.n_blocks_dim];
        status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
        mkl_sparse_status(status);
    }

    return;
}