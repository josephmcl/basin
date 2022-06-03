#pragma once

#include <cmath>
#include <numbers>
#include <vector>
#include <tuple>
#include <iostream>
#include <array>

#include "definitions.h"
#include "ranges.h"
#include "numerical.h"
#include "solve.h"


/* petsc headers */
#include "petscsys.h"
#include "petscmat.h" 
#include "petscviewer.h"
#include "petscvec.h"
#include "petscksp.h"

namespace sbp_sat { 


  using nat_t = std::size_t;
  using real_t = type::real_t;  
  using vector_t = std::vector<real_t>;
  using real_v = std::vector<real_t>;
  using range_t = numerics::linrange<real_t>;

  using block_t = std::tuple<range_t, range_t>; 
  using block_v = std::vector<block_t>;

  using domain_t = std::tuple<real_t, real_t>;
  using domain_v = std::vector<domain_t>;

  using boundary_t = std::tuple<range_t, nat_t>;
  using boundary_tx2 = std::array<std::tuple<range_t, nat_t>, 4>;
  using boundary_vx2 = std::vector<boundary_tx2>;

  auto const π = std::numbers::pi_v<real_t>;
  auto static to_real_t = [](nat_t n){return static_cast<real_t>(n);};
  using petsc_matrix = Mat;
  using petsc_vector = Vec;


namespace x2 {


    using ℤ = std::size_t;
    using ℝ = type::real_t;  

    /* */
    template <typename T>
    T analytical_solution(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return std::sin(cx * pi * x + cy * pi * y); 
    }

    template <typename T>
    T uxx(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return -(cx * cx) * (pi * pi) * std::sin(cx * pi * x + cy * pi * y); 
    }

    template <typename T>
    T uyy(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return -(cy * cy) * (pi * pi) * std::sin(cx * pi * x + cy * pi * y); 
    }

    template <typename T>
    T ux(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return cx * pi * std::cos(cx * pi * x + cy * pi * y); 
    }

    template <typename T>
    T uy(T &x, T &y, T const &cx=1., T const &cy=1.) {
        auto constexpr pi = std::numbers::pi_v<T>;
        return cy * pi * std::cos(cx * pi * x + cy * pi * y); 
    }

    void petsc_hybridized_poisson(
      real_v             &result,
      domain_v     const &domain, 
      block_v      const &blocks, 
      boundary_vx2 const &boundary);

    void write_m( 
      petsc_matrix       &m,
      block_t      const &block_x1, 
      block_t      const &block_x2, 
      boundary_tx2 const &boundary_x1, 
      boundary_tx2 const &boundary_x2);

    void write_d2_h1(
      petsc_matrix                   &M, 
      std::vector<long double> const &h, 
      range_t                  const &local,
      long double              const  spacing_square); 

}; /* namespace sbp_sat::x2 */
}; /* namespace sbp_sat     */
