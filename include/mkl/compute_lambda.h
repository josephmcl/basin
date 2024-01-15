
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
    long long int *piv,
    components &sbp);

void compute_lambda(
    real_t *lambdaA,
    long long int *piv,
    real_t *lambdab,
    components &sbp);