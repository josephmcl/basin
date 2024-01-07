#pragma once 

#include <cmath>
#include <numbers>
#include <vector>
#include <tuple>
#include <iostream>
#include <array>
#include <optional>
#include <functional>

#include "definitions.h"
#include "numerical.h"
#include "csr.h"

/* petsc headers */

namespace sbp_sat_2d { 

  using real_t = type::real_t; 
  using vector_t = std::vector<real_t>;
  auto const π = std::numbers::pi_v<real_t>;
  auto static to_real_t = [](std::size_t n){return static_cast<real_t>(n);};

  // Make and store all the components needed to do SBP-SAT things. 
  class components {

  public:
    components(
      std::size_t const points,
      real_t const span,
      std::size_t const accuracy = 2);
    ~components();

    std::size_t const n;
    real_t      const span;
    std::size_t const accuracy;
    real_t      const spacing;

    real_t τ, β;

    std::size_t n_blocks, n_interfaces;
    std::size_t n_blocks_dim; 

    std::size_t n_threads;

    /*
    petsc_matrix hl;
    petsc_matrix hx, hy;
    petsc_matrix h1x, h1y;
    petsc_matrix bsx, bsy;  
    petsc_matrix ln, ls, le, lw;  
    petsc_matrix d2x, d2y;

    std::vector<real_t> h1v;

  private: 
    void make_bs();
    void make_l();
    void make_h();
    void make_h1();
    void make_hl();
    void make_d2();
    */
	};
}; // namespace sbp_sat_2d