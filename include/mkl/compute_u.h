#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_sparse_qr.h"


#include "omp.h"

void compute_u(
    real_t *u,
   std::vector<double *> &M,
   std::vector<int *> &Mpiv,
    real_t *rhs, 
    components &sbp);