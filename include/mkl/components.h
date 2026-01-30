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


using real_t = type::real_t; 
using vector_t = std::vector<real_t>;
auto const PI_VALUE = 3.14159;
auto static to_real_t = [](std::size_t n){return static_cast<real_t>(n);};

// Make and store all the components needed to do SBP-SAT things. 
class components {

public:
  components(
    std::size_t const points,
    real_t const span,
    std::size_t const accuracy = 2);

  std::size_t const n;
  real_t      const span;
  std::size_t const accuracy;
  real_t      const spacing;

  real_t TAU_VALUE, BETA_VALUE;

  std::size_t n_blocks, n_interfaces;
  std::size_t n_blocks_dim; 

  std::size_t n_threads;
  std::size_t n_ranks;
  std::size_t rank;

  std::size_t rank_limit_u;
  std::vector<std::size_t> rank_index_u;

  std::size_t rank_limit_iu;
  std::vector<std::size_t> rank_index_iu;
  
  csr<real_t> hl;
  csr<real_t> hx, hy;
  csr<real_t> h1x, h1y;
  csr<real_t> bsx, bsy;  
  csr<real_t> ln, ls, le, lw;  
  csr<real_t> d2x, d2y;

  std::vector<real_t> h1v;

private: 
  void make_bs();
  void make_l();
  void make_h();
  void make_h1();
  void make_hl();
  void make_d2();
};