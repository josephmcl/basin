#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_sparse_qr.h"

void compute_u(
    real_t *u,
    std::vector<sparse_matrix_t> &M,
    real_t *rhs, 
    components &sbp);