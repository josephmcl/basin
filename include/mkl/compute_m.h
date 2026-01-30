#pragma once 

#include <array>

#include "components.h"
#include "error.h"

#include "mkl.h"
#include "mkl_spblas.h"


void make_m(
  sparse_matrix_t                  *M,
  components                       &sbp, 
  std::array<std::size_t, 4> const &boundary);

void make_m_boundary(
  sparse_matrix_t *M,
  components const &sbp, 
  std::size_t const direction,
  std::size_t const boundary);