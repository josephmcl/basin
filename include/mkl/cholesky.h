#pragma once

#define from(a,b) std::make_tuple(a,b)

/* STL headers       */
#include <tuple>
#include <optional>

/* Project headers   */
#include "components.h"
#include "csr.h"
#include "error.h"

/* MKL headers       */
#include "mkl.h"

// TODO: move this to definitions.h
using error = std::optional<int>;

namespace cholesky {

    /* Wrapper for LAPACK Cholesky factorize.
       Overwrites *a. */
    error factor(
        double     *a, 
        components &sbp,
        size_t      n = 0);

    /* Wrapper for LAPACK Cholesky factorize.
       Overwrites *a. */
    error factor_rfp(
        double     *a, 
        components &sbp);
        
    /* Wrapper for LAPACK Cholesky solver.
       Overwrites *b. */
    /*error solve(
        double     *a, 
        double     *b, 
        components &sbp);*/

    /*error solve_rfp(
        double     *a, 
        double     *b, 
        components &sbp);*/

    inline error solve(
        double     *a, 
        double     *b, 
        components &sbp,
        std::size_t n = 0) {

        constexpr int layout = LAPACK_COL_MAJOR;
        constexpr char ul = 'L';
        constexpr int nrhs = 1;
        if (n == 0) n = sbp.n * sbp.n;
        lapack_int res = LAPACKE_dpotrs(
            layout, ul, n, nrhs, a, n, b, n);
        return (res == 0) 
            ? std::nullopt : error(res);
    }

    inline error solve_rfp(
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

    // NOTE: const isn't used pretty much 
    //       anywhere because it causes 
    //       issues with MKL. 

} /* end namespace cholesky */

using sparse_matrix_z = sparse_matrix_t;
using dense_matrix_z = double *;

void to_dense (
    dense_matrix_z    &b, 
    sparse_matrix_z   &a,
    std::size_t        n);

void to_rfp(
    dense_matrix_z    &b, 
    sparse_matrix_z   &a,
    std::size_t        n);
