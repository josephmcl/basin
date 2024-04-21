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

void compute_mg(
    real_t *Mg,
    std::vector<double *> &M,
    std::vector<int *> &Mpiv,
    real_t *g,
    components &sbp);

void compute_mg_mpi(
    real_t *Mg,
    vv<sparse_matrix_t> &M,
    real_t *g,
    components &sbp);
