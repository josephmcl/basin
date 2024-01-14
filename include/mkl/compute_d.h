#pragma once 

#include <vector>

#include "components.h"
#include "error.h"
#include "csr.h"
#include "definitions.h"

#include "mkl.h"
#include "mkl_spblas.h"

using real_t = type::real_t;

void compute_d(
  sparse_matrix_t          *D, 
  components               &sbp, 
  vv<std::size_t> const &interfaces);