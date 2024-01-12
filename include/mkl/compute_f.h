#pragma once 

#include <vector>

#include "components.h"
#include "error.h"
#include "csr.h"
#include "definitions.h"

#include "mkl.h"
#include "mkl_spblas.h"

using real_t = type::real_t;

void compute_f(
  std::vector<sparse_matrix_t> &F_sparse,
  std::vector<real_t *>        &F_dense,
  components                   &sbp);


void fcompop(
  sparse_matrix_t *f_sparse, 
  real_t *f_dense,
  csr<real_t> &L, 
  csr<real_t> &B,
  csr<real_t> &H,
  real_t       const τ, 
  real_t       const β);
  