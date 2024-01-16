
#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_sparse_qr.h"

void initialize_lambda(
    double *lambdaA,
    MKL_INT *piv,
    components &sbp);

void compute_lambda(
    real_t *lambdaA,
    MKL_INT *piv,
    real_t *lambdab,
    components &sbp);