#pragma once

#include "error.h"
#include "components.h"
#include "csr.h"
#include "definitions.h"
#include "error.h"
#include "timing.h"

#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_sparse_qr.h"

#include <cstring>

#include "omp.h"

using real_t = type::real_t;

// std::size_t compute_mf(
//  std::vector<real_t *>        &x,
//  vv<sparse_matrix_t> &m,
//  std::vector<real_t *>        &f,
//  components             const &sbp);

using mf_instrument = timing::instrument< 
  std::vector<real_t *>  &,
  vv<sparse_matrix_t>    &,
  std::vector<real_t *>  &,
  components      const  &>;

extern mf_instrument compute_mf;
