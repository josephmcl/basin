#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include "mkl.h"
#include "mkl_spblas.h"

void compute_lambda_a(
  sparse_matrix_t              *lambdaA,
  sparse_matrix_t              *D,
  std::vector<sparse_matrix_t> &F,
  std::vector<real_t *>        &MF,
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components                   &sbp);

