#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include "mkl.h"
#include "mkl_spblas.h"

void compute_rhs(
    real_t *rhs,
    std::vector<sparse_matrix_t> &F,
    real_t *lambda,
    vv<std::size_t> &F_symbols,
    components &sbp);