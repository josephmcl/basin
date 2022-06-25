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
#include "ranges.h"
#include "numerical.h"
#include "solve.h"

#include "linalg.h"


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

  using namespace linalg;

  auto constexpr fw = framework::petsc;

  template <typename T>
  using optref = std::optional<std::reference_wrapper<T>>;

  auto constexpr transposed = true;

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
    petsc_matrix       &M, 
    real_v       const &h, 
    range_t      const &local,
    real_t       const spacing_square, 
    real_t       const coeff=1.); 

  void write_h1(
    petsc_matrix       &M, 
    real_v       const &h); 

void write_boundary_x1_left(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bsx1,
  vector_t     const &hx2, 
  long double  const  β,
  long double  const  τ,
  int          const  order = 1);

void write_boundary_x1_right(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bsx1,
  vector_t     const &hx2, 
  real_t       const  β,
  real_t       const  τ,
  int          const  order = 1); 

void write_boundary_x2_left(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bsx2,
  vector_t     const &hx1, 
  real_t       const  β,
  real_t       const  τ,
  int          const  order = 1);

void write_boundary_x2_right(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bsx2,
  vector_t     const &hx1, 
  real_t       const  β,
  real_t       const  τ,
  int          const  order = 1);


/* 

            (1, 0)           
            ______
           |      |
    (0, 0) |      | (0, 1)
           |______|

            (1, 1)

   _______      _______      _______      _______
  |\      |    |       |    |`      |    |,      |
  |  \    |    |       |    |  `    |    |  ,    |
  |       |    |    \  |    |    `  |    |    ,  |
  |_______|    |______\|    |______`|    |______,|
    
    (0,0)        (0,1)        (1,0)        (1,1)

*/

/* Set a partial Kronecker product i.e., lower or upper, for a single 
   diagonal matrix. Where the full Kronecker product would result in 
   entries on a single boundary row or column. */
template<nat_t major, nat_t minor> 
void boundary_diagonal(
  matrix<fw>       &dest,
  nat_t      const  size1,
  nat_t      const  size2,
  vector_t   const &h, 
  real_t     const  c = 1.) {

  nat_t row, col; 
  real_t val;
  nat_t const size = size1 * size2;
  if constexpr (major == 0 and minor == 0) {
    for (std::size_t i = 0; i < h.size(); ++i) {
      row = i;
      col = i;
      val = c * h[i];
      set_matrix_value<fw>(dest, row, col, val);
    }
  }
  else if constexpr (major == 0 and minor > 0) {
    for (std::size_t i = 0; i < h.size(); ++i) {
      row = size - i - 1;
      col = size - i - 1;
      val = c * h[i];
      set_matrix_value<fw>(dest, row, col, val);
    } 
  }
  else if constexpr (major > 0 and minor == 0) {
    for (std::size_t i = 0; i < h.size(); ++i) {
      col = i * size2;
      val = c * h[i];
      row = i * size2;
      set_matrix_value<fw>(dest, row, col, val);
    } 
  }
  else if constexpr (major > 0 and minor > 0) {
    for (std::size_t i = 0; i < h.size(); ++i) { 
      col = size - (i * size2);
      val = c * h[i];
      row = size - (i * size2);
      set_matrix_value<fw>(dest, row, col, val);
    } 
  }
}


template<std::size_t major, std::size_t minor, bool transpose = false> 
void boundary_bs(
  matrix<fw>             & dest,
  nat_t            const   size1,
  nat_t            const   size2,
  vector_t         const & bs, 
  optref<vector_t  const>  h = {},
  real_t           const   c = 1.) {

  std::function<real_t(nat_t, nat_t)> entry;
  if (h) {
    vector_t tmp = *h;
    entry = [&bs, h = tmp, &c] (nat_t i, nat_t j) -> real_t {
      return c * h[i] * bs[j];};
  }
  else
    entry = [&bs, &c] (nat_t i, nat_t j) -> real_t {
        (void) i; return c * bs[j];};

  nat_t row, col; 
  real_t val;
  nat_t const size = size1 * size2;
  if constexpr (major == 0 and minor == 0 and transpose) {
    for (std::size_t i = 0; i < size2; ++i) {
      col = i;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        row = i + (j * size1);
        set_matrix_value<fw>(dest, row, col, val);
      }
    }
  }
  else if constexpr (major == 0 and minor == 0 and not transpose) {
    for (std::size_t i = 0; i < size2; ++i) {
      row = i;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        col = i + (j * size1);
        set_matrix_value<fw>(dest, row, col, val);
      }
    }
  }
  else if constexpr (major == 0 and minor > 0 and transpose) {
    for (std::size_t i = 0; i < size2; ++i) {
      col = size - i;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        row = size - i - (j * size1);
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major == 0 and minor > 0 and not transpose) {
    for (std::size_t i = 0; i < size2; ++i) {
      row = size - i;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        col = size - i - (j * size1);
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major > 0 and minor == 0 and transpose) {
    for (std::size_t i = 0; i < size1; ++i) {
      col = i * size2;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        row = i + j;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major > 0 and minor == 0 and not transpose) {
    for (std::size_t i = 0; i < size1; ++i) {
      row = i * size2;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        col = i + j;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major > 0 and minor > 0 and transpose) {
    for (std::size_t i = 0; i < size1; ++i) { 
      col = size - (i * size2);
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        row = size - i - j;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major > 0 and minor > 0 and not transpose) {
    for (std::size_t i = 0; i < size1; ++i) { 
      row = size - (i * size2);
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        col = size - i - j;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
}

}; /* namespace sbp_sat::x2 */
}; /* namespace sbp_sat     */
