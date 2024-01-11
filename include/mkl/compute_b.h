#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"

#include "mkl.h"
#include "mkl_spblas.h"

using real_t = type::real_t;

void 
compute_b(
    std::vector<sparse_matrix_t>    &  B, 
    components                      &sbp);

void compute_b1(
    sparse_matrix_t    & B, 
    csr<real_t>        & H, 
    real_t       const   τ, 
    csr<real_t>        & L, 
    real_t       const   β, 
    csr<real_t>        &BS,
    std::size_t  const   n);

void compute_b2(
    sparse_matrix_t    & B, 
    csr<real_t>        & H, 
    real_t       const   τ, 
    csr<real_t>        & L, 
    real_t       const   β, 
    csr<real_t>        &BS,
    std::size_t  const   n);