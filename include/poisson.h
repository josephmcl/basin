#pragma once

#include <cmath>
#include <numbers>
#include <vector>
#include <tuple>
#include <iostream>

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

/* setup the Poisson Equation in 2D */
namespace poisson1d {

  using nat_t = std::size_t;
  using real_t = type::real_t;  
  using vector_t = std::vector<real_t>;
  using range_t = numerics::linrange<real_t>;
  using domain_t = std::tuple<real_t, real_t>;
  using boundary_t = std::tuple<real_t, nat_t, real_t, nat_t>;
  auto const π = std::numbers::pi_v<real_t>;
  auto static to_real_t = [](nat_t n){return static_cast<real_t>(n);};
  using petsc_matrix = Mat;
  using petsc_vector = Vec;

  void analytical(
    vector_t       &res,
    nat_t  const  size, 
    real_t const  left, 
    real_t const  right);

  void analytical_blocked(
    vector_t       &res,
    nat_t    const  size, 
    nat_t    const  local_problems,
    real_t   const  left = 0, 
    real_t   const  right= 1);

  domain_t   const domain   = std::make_tuple(0., 1.);
  boundary_t const boundary = std::make_tuple(0., 1, -π, 2);

  void petsc_problem(
    vector_t         &result,
    nat_t      const  size, 
    nat_t      const  local_problems = 2,
    domain_t   const  domain         = domain, 
    boundary_t const  boundary       = boundary);

  void petsc_hybridized_problem(
    vector_t         &result,
    nat_t      const  size, 
    nat_t      const  local_problems = 2,
    domain_t   const  domain         = domain, 
    boundary_t const  boundary       = boundary);

  /*  Set the second-order SBP D2 operator in A using D2 sized to the 
      local_domain_size and partitioned by the number of 
      local_domain_ranges. */
  void write_d2(
    petsc_matrix               &A, 
    std::vector<range_t> const &local_domain_ranges,
    long double          const  local_domain_size,
    long double          const  spacing_square);

  void write_d2_h1(
    petsc_matrix                   &M, 
    std::vector<long double> const &h1, 
    std::vector<range_t>     const &local_domain_ranges,
    long double              const  local_domain_size,
    long double              const  spacing_square); 

  void write_boundaries(
    petsc_matrix               &A, 
    std::vector<range_t> const &local_domain_ranges,
    vector_t             const &hi, 
    vector_t             const &bs, 
    long double          const  matrix_points,
    long double          const  local_domain_size,
    long double          const  β,
    long double          const  σ1,
    long double          const  σ2);

  void write_fluxs(
    petsc_matrix               &A, 
    std::vector<range_t> const &local_domain_ranges,
    vector_t             const &hi, 
    vector_t             const &bs, 
    long double          const  local_domain_size,
    long double          const  β,
    long double          const  σ1,
    long double          const  ϵ);

  void write_upper_first_order_boundary(
    petsc_matrix               &M, 
    vector_t             const &bs, 
    long double          const  β,
    long double          const  τ);

  void write_lower_first_order_boundary(
    petsc_matrix               &M, 
    vector_t             const &bs,
    long double          const  n, 
    long double          const  β,
    long double          const  τ);

  void write_lower_second_order_boundary(
    petsc_matrix               &M, 
    vector_t             const &bs,
    long double          const  n, 
    long double          const  β,
    long double          const  τ);
};
