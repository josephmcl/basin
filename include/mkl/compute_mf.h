#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"

#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_pardiso.h"

#include <cstring>

using real_t = type::real_t;

void compute_mf(
  std::vector<real_t *>        &x,
  std::vector<sparse_matrix_t> &m,
  std::vector<real_t *>        &f,
  components             const &sbp);