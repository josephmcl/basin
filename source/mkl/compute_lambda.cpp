
#include "compute_lambda.h"
#include "omp.h"

void initialize_lambda(
    double *a,
    MKL_INT *piv,
    components &sbp) {
        
    
    MKL_INT size = static_cast<MKL_INT>(sbp.n * sbp.n_interfaces);
    //std::cout << size << std::endl;
    //std::cout << size * size * sizeof(double) << std::endl;
    MKL_INT err;
    
    /*a
    // double *aa = (double *) malloc(size * size * sizeof(double));
    // memset(aa, 1, size * sizeof(double));

    std::cout << "tomato" << std::endl;

    
        

        for (std::size_t j = 0; j != sbp.n_interfaces * sbp.n * sbp.n_interfaces * sbp.n; ++j) {
            aa[j] = 1.;
        }

        // std::size_t v;
        for (std::size_t j = 0; j != sbp.n_interfaces * sbp.n; ++j) {
            //v = j % sbp.n;
            // lambdaA[(j * sbp.n_interfaces * sbp.n) + j] += sbp.h1v[v] * 2. * sbp.τ;
            std::cout << (j * sbp.n_interfaces * sbp.n) + j << std::endl;
            aa[(j * sbp.n_interfaces * sbp.n) + j] = 2;
        }
    */    
    //#pragma omp parallel num_threads(1)
    //{
    // char u = 'U';
    //dpotrf2(&u, &size, a, &size, &err);

    // auto tds = omp_get_num_threads();
    // omp_set_num_threads(56);
    // mkl_set_num_threads(56);

    dgetrf(&size, &size, a, &size, piv, &err);
    if (err != 0)
        std::cout << "dgetrf err: " << err << std::endl;

    // omp_set_num_threads(tds);
    // mkl_set_num_threads(tds);
}

void compute_lambda(
    real_t *lambdaA,
    MKL_INT *piv,
    real_t *lambda,
    components &sbp) {

    //auto tds = omp_get_num_threads();
    //mkl_set_num_threads(56);

    char t = 'N';
    MKL_INT size = static_cast<MKL_INT>(sbp.n * sbp.n_interfaces);
    //std::cout << size << std::endl;
    MKL_INT rhs = 1;
    MKL_INT s;
    dgetrs(
        &t, &size, 
        &rhs, lambdaA, &size,  piv, 
        lambda, &size, &s);
    if (s != 0)
            std::cout << "dgetrs err: " << s << std::endl;

    //mkl_set_num_threads(tds);

    /*
    sparse_status_t status;
    status = mkl_sparse_d_qr_solve(
        SPARSE_OPERATION_NON_TRANSPOSE, *lambdaA, nullptr,
        SPARSE_LAYOUT_COLUMN_MAJOR, 1, lambda, sbp.n * sbp.n_interfaces, 
        lambdab, sbp.n * sbp.n_interfaces);
    mkl_sparse_status(status);
    */
}
