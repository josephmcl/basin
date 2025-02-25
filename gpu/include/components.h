#pragma once 

#include <cmath>
//#include <numbers>
#include <vector>
#include <tuple>
#include <iostream>
#include <array>
#include <optional>
#include <functional>

#include "definitions.h"
#include "numerical.h"
#include "csr.h"

#define using_rocblas 1
#ifdef using_rocblas
  #include "rocblas.h"
  using rocblas_handle_type = rocblas_handle;
#else 
  using rocblas_handle_type = void *;
#endif

// Make and store all the components needed to do SBP-SAT things. 
class components {

public:
  components(
    std::size_t const points,
    real_t const span,
    std::size_t const accuracy = 2);

  std::size_t const n = 0;
  real_t      const span = 0;
  std::size_t const accuracy = 0;
  real_t      const spacing = 0;

  real_t 𝜏, β;

  std::size_t n_blocks, n_interfaces;
  std::size_t n_blocks_dim; 

  std::size_t n_threads;
  std::size_t n_ranks;
  std::size_t rank;

  std::size_t rank_limit_u;
  std::vector<std::size_t> rank_index_u;

  
  csr<real_t> hl;
  csr<real_t> hx, hy;
  csr<real_t> h1x, h1y;
  csr<real_t> bsx, bsy;  
  csr<real_t> ln, ls, le, lw;  
  csr<real_t> d2x, d2y;

  std::vector<real_t> h1v;

  rocblas_handle_type rb_handle;

private: 
  void make_bs();
  void make_l();
  void make_h();
  void make_h1();
  void make_hl();
  void make_d2();
};