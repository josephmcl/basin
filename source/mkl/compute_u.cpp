#include "compute_u.h"

//#include "ittnotify.h"
#include "timing.h"
#include "stdio.h"

/*
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

    #pragma omp parallel for private(up, rp, k)
    for (std::size_t i = 0; i < sbp.rank_limit_u; ++i) {
       
        auto td = omp_get_thread_num();
        auto ii = sbp.rank_index_u[i];
        up = &u[i * sbp.n * sbp.n];
        rp = &rhs[i * sbp.n * sbp.n];
        // rank_index_u is the global index, thats what we modulo. 
        k = mi[ii % sbp.n_blocks_dim];
        
        // mkl_set_num_threads_local(4);
        status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, up , sbp.n * sbp.n, 
            rp, sbp.n * sbp.n);
        mkl_sparse_status(status);
    }    
    return;
}
*/


void compute_u_sqr(
    real_t *Mg,
    vv<sparse_matrix_t> &M,
    real_t *g,
    components &sbp) {  

    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

    //__itt_domain* domain = __itt_domain_create("Compute U Domain");
    //__itt_string_handle* task = __itt_string_handle_create("Compute U Task");

    sparse_status_t status;
    double *gp, *mgp;
    std::size_t k;
    std::size_t limit = sbp.rank_limit_u;   

    //__itt_task_begin(domain, __itt_null, __itt_null, task);

    // auto begin = timing::read();

    
    #pragma omp parallel for private(gp, mgp, k) 
    for (std::size_t i = 0; i < limit; ++i) {
        auto td = omp_get_thread_num();
        auto ii = sbp.rank_index_u[i];
        gp = &g[i * sbp.n * sbp.n];
        mgp = &Mg[i * sbp.n * sbp.n];
        // rank_index_u is the global index, thats what we modulo. 
        k = mi[ii % sbp.n_blocks_dim];

        //std::cout << td << " " << std::endl;

        // status = 
        // auto begin = timing::read();
        //mkl_set_num_threads_local(8);
        mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
        // auto end = timing::read();
        // ::printf("%f\n", end - begin);
        // mkl_sparse_status(status);
    }

    //auto end = timing::read();
    //std::cout << "close timer " << end - begin << std::endl;

    //__itt_task_end(domain);
}

void compute_u_dc(
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

void compute_u_rfpc(
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