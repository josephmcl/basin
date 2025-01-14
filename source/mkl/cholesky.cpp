#include "cholesky.h"


error cholesky::factor(
        double     *a, 
        components &sbp,
        std::size_t n) {
    constexpr int layout = LAPACK_COL_MAJOR;
    constexpr char ul = 'L';
    if (n == 0)
        n = sbp.n * sbp.n;
    lapack_int res = LAPACKE_dpotrf(
        layout, ul, n, a, n);
    return (res == 0) 
        ? std::nullopt : error(res);
}

error cholesky::factor_rfp(
        double     *a, 
        components &sbp) {
    constexpr int layout = LAPACK_COL_MAJOR;
    constexpr char transr = 'N';
    constexpr char ul = 'L';
    const std::size_t n = sbp.n * sbp.n;
    lapack_int res = LAPACKE_dpftrf(
        layout, transr, ul, n, a);
    return (res == 0) 
        ? std::nullopt : error(res);
}

/*

error cholesky::solve(
        double     *a, 
        double     *b, 
        components &sbp) {
    constexpr int layout = LAPACK_COL_MAJOR;
    constexpr char ul = 'L';
    constexpr int nrhs = 1;
    const std::size_t n = sbp.n * sbp.n;
    lapack_int res = LAPACKE_dpotrs(
        layout, ul, n, nrhs, a, n, b, n);
    return std::nullopt;
}


inline error cholesky::solve_rfp(
        double     *a, 
        double     *b, 
        components &sbp) {
    constexpr int layout = LAPACK_COL_MAJOR;
    constexpr char transr = 'N';
    constexpr char ul = 'L';
    constexpr int nrhs = 1;
    const std::size_t n = sbp.n * sbp.n;
    // mkl_set_num_threads_local(4);
    lapack_int res = LAPACKE_dpftrs(
        layout, transr, ul, n, nrhs, a, b, n);
    return (res == 0) 
        ? std::nullopt : error(res);
    // return std::nullopt; // for prod
}
*/


void to_dense (
        dense_matrix_z  &b,
        sparse_matrix_z &a,
        std::size_t      n) {
            
    std::size_t const size 
        = sizeof(double) * n * n;
    b = (double *) mkl_malloc(size, 64);

    csr<double> eye(n, n);
    for (std::size_t i = 0; i != n; ++i) {
        eye(1., i, i);
    }
    sparse_matrix_t i;
    eye.mkl(&i);

    sparse_operation_t constexpr transpose 
        = SPARSE_OPERATION_NON_TRANSPOSE;
    sparse_layout_t constexpr layout 
        = SPARSE_LAYOUT_ROW_MAJOR;
    sparse_status_t status =  mkl_sparse_d_spmmd(
        transpose, a, i, layout, b, n);
    mkl_sparse_status(status);

    return;
}

void to_rfp (
        dense_matrix_z  &b,
        sparse_matrix_z &a,
        std::size_t      n) {
            
    std::size_t size 
        = sizeof(double) * n * n;
    double *temp = (double *)
        mkl_malloc(size, 64);
    

    /* Make an indentity matrix. */
    csr<double> eye(n, n);
    for (std::size_t i = 0; i != n; ++i)
        eye(1., i, i);
    sparse_matrix_t i;
    eye.mkl(&i);

    /* Make a full dense matrix. */
    sparse_operation_t constexpr transpose 
        = SPARSE_OPERATION_NON_TRANSPOSE;
    sparse_layout_t constexpr layout 
        = SPARSE_LAYOUT_COLUMN_MAJOR;
    sparse_status_t status =  mkl_sparse_d_spmmd(
        transpose, a, i, layout, temp, n);
    mkl_sparse_status(status);

    /* Allocate storage for rfp format. */
    size = sizeof(double) * n * (n + 1) / 2;
    b = (double *) mkl_malloc(size, 64);

    char constexpr transr = 'N';
    char constexpr uplo = 'L';
    lapack_int res = LAPACKE_dtrttf(
        layout, transr, uplo, n, temp, n, b);

    mkl_free(temp);

    if (res != 0) {
        std::cerr << "LAPACK routine dtrttf "
            << "returned error " << res
            << std::endl;
    }

    return;
}