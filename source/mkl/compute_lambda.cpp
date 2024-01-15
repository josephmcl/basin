
#include "compute_lambda.h"

void initialize_lambda(
    double *a,
    long long int *piv,
    components &sbp) {
    
    std::size_t size = sbp.n * sbp.n_interfaces;
    lapack_int s = LAPACKE_dgetrf(
        LAPACK_ROW_MAJOR, size, 
        size, a, size, piv);
    std::cout << s << std::endl;
    /*
    

    /
    matrix_descr da;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    sparse_status_t status;
    status = mkl_sparse_qr_reorder(*lambdaA, da);
    mkl_sparse_status(status);
    std::cout << "reordered" << std::endl;
    status = mkl_sparse_d_qr_factorize(*lambdaA, nullptr);
    mkl_sparse_status(status);
    */
}

void compute_lambda(
    real_t *lambdaA,
    long long int *piv,
    real_t *lambdab,
    components &sbp) {

    lapack_int s = LAPACKE_dgetrs(
        LAPACK_ROW_MAJOR, 'N', sbp.n * sbp.n_interfaces, 
        1, lambdaA, sbp.n * sbp.n_interfaces, piv, 
        lambdab, sbp.n * sbp.n_interfaces);
    std::cout << s << std::endl;
    /*
    sparse_status_t status;
    status = mkl_sparse_d_qr_solve(
        SPARSE_OPERATION_NON_TRANSPOSE, *lambdaA, nullptr,
        SPARSE_LAYOUT_COLUMN_MAJOR, 1, lambda, sbp.n * sbp.n_interfaces, 
        lambdab, sbp.n * sbp.n_interfaces);
    mkl_sparse_status(status);
    */
}
