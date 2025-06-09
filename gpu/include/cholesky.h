#pragma once

#define from(a,b) std::make_tuple(a,b)

/* STL headers       */
#include <tuple>
#include <optional>
#include <climits>

/* Project headers   */
#include "components.h"
#include "csr.h"
#include "error.h"

#include "rocblas.h"
#include "rocsolver.h"


// TODO: move this to definitions.h
using error = std::optional<int>;


namespace cholesky {
        
    /* Wrapper for LAPACK Cholesky solver.
       Overwrites *b. */

    inline error solve(
        double     *a, 
        double     *b, 
        components &sbp,
        std::size_t n = 0) {

        constexpr int nrhs = 1;
        constexpr rocblas_fill uplo = rocblas_fill_lower;
        if (n == 0) n = sbp.n * sbp.n;
        rocblas_status res = rocsolver_dpotrs(
            sbp.rb_handle, uplo, n, nrhs, a, n, b, n);
        hipDeviceSynchronize();
        return (res == 0) ? std::nullopt : error(res);
    } 

    inline error solve_no_sync(
        double     *a, 
        double     *b, 
        components &sbp,
        std::size_t n = 0) {

        constexpr int nrhs = 1;
        constexpr rocblas_fill uplo = rocblas_fill_lower;
        if (n == 0) n = sbp.n * sbp.n;
        rocblas_status res = rocsolver_dpotrs(
            sbp.rb_handle, uplo, n, nrhs, a, n, b, n);
        return (res == 0) ? std::nullopt : error(res);
    } 

} /* end namespace cholesky */
