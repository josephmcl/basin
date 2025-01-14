#include "compute_mg.h"

//#include "ittnotify.h"


void compute_mg_sqr(
    real_t *Mg,
    vv<sparse_matrix_t> &M,
    real_t *g,
    components &sbp) {  

    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

    sparse_status_t status;
    double *gp, *mgp;
    std::size_t k;
    std::size_t limit = sbp.rank_limit_u;   

    //__itt_domain* domain = __itt_domain_create("Compute MG Domain");
    //__itt_string_handle* task = __itt_string_handle_create("Compute MG Task");

    std::cout << omp_get_thread_num() << std::endl; 


    //__itt_task_begin(domain, __itt_null, __itt_null, task);

    // NOTE: This code needs to be called from within an OpenMP 
    //       parallel region. 
    #pragma omp for private(gp, mgp, k)
    for (std::size_t i = 0; i < limit; ++i) {
        auto td = omp_get_thread_num();
        auto ii = sbp.rank_index_u[i];
        gp = &g[i * sbp.n * sbp.n];
        mgp = &Mg[i * sbp.n * sbp.n];
        // rank_index_u is the global index, thats what we modulo. 
        k = mi[ii % sbp.n_blocks_dim];

        status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
        mkl_sparse_status(status);
    }

    //__itt_task_end(domain);
}

void compute_mg_dc(
    real_t                *Mg,
    std::vector<double *> &M,
    components            &sbp) {  

    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

    sparse_status_t status;
    double *mgp;
    std::size_t k;
    std::size_t limit = sbp.rank_limit_u;   
    #pragma omp parallel for private(mgp, k) 
    for (std::size_t i = 0; i < limit; ++i) {
        auto td = omp_get_thread_num();
        auto ii = sbp.rank_index_u[i];
        mgp = &Mg[i * sbp.n * sbp.n];
        // rank_index_u is the global index, thats what we modulo. 
        k = mi[ii % sbp.n_blocks_dim];

        auto error = cholesky::solve(M[k], mgp, sbp);
        if (error) {
            std::cout << "Choleksy solve failed code " 
                << *error << std::endl;
        }
        
    }
}

void compute_mg_rfpc(
    real_t                *Mg,
    std::vector<double *> &M,
    components            &sbp) {  

    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

    sparse_status_t status;
    double *mgp;
    std::size_t k;
    std::size_t limit = sbp.rank_limit_u;   
    #pragma omp parallel for private(mgp, k) 
    for (std::size_t i = 0; i < limit; ++i) {
        auto td = omp_get_thread_num();
        auto ii = sbp.rank_index_u[i];
        mgp = &Mg[i * sbp.n * sbp.n];
        // rank_index_u is the global index, thats what we modulo. 
        k = mi[ii % sbp.n_blocks_dim];

        auto error = cholesky::solve_rfp(M[k], mgp, sbp);
        
        if (error) {
            std::cout << "Choleksy solve failed code " 
                << *error << std::endl;
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
    #pragma omp parallel for private(gp, mgp, k)
    for (std::size_t i = 0; i < limit; ++i) {
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

