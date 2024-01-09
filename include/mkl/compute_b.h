#pragma once

#include "components.h"
#include "csr.h"
#include "definitions.h"

#include "mkl.h"
#include "mkl_spblas.h"

using real_t = type::real_t;

void 
compute_b(
    std::vector<csr<real_t>>        &  B, 
    components                      &sbp);

void compute_b1(
    csr<real_t>        & B, 
    csr<real_t>        & H, 
    real_t       const   τ, 
    csr<real_t>        & L, 
    real_t       const   β, 
    csr<real_t>        &BS,
    std::size_t  const   n);

void compute_b2(
    csr<real_t>        & B, 
    csr<real_t>        & H, 
    real_t       const   τ, 
    csr<real_t>        & L, 
    real_t       const   β, 
    csr<real_t>        &BS,
    std::size_t  const   n);