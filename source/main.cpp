#include <ranges>
#include <vector>
#include <algorithm>
#include <tgmath.h>
#include <iomanip>
#include <optional>
#include <iostream>
#include <cmath>

#include "infrastructure.h" 
#include "ranges.h"
#include "domain.h"
#include "physical.h"
#include "numerical.h"
#include "definitions.h"
#include "poisson2d.h"

#include "unistd.h"

// #include "plot.h"

#include "solve.h"
#include "poisson.h"
#include "hybrid_sbp_sat_2d.h"

#include "timing.h"

void poisson_1d_convergence_test(
  std::size_t volume_points, 
  std::size_t blocks);

void poisson_1d_hybrid_convergence_test(
  std::size_t volume_points, 
  std::size_t blocks);

void main_2d();

int main (int argc, char *argv[]) {

  // Initialize petsc, among other libraries.
  if (!infrastructure::initialize()) exit(-1);; { 

    // timing::init();
    // poisson_1d_hybrid_convergence_test(199, 3);

    main_2d();

  // Cleanup petsc 
  } infrastructure::cleanup(); return 0;
}

void poisson_1d_convergence_test(
  std::size_t volume_points, 
  std::size_t blocks) {

  /* Solving the Poisson Equation in 1D

    ⎧ u_xx π^2 sin(πx) = 0 , 0 < x < 1
    ⎨ u = 0    (Dirichlet) , x = 0 
    ⎩ u_x = -π (Neumann)   , x = 1 

  */

  std::cout << "[ Poisson 1D equation with arbitrary flux " 
    << "interfaces. ]" << std::endl;

  std::cout << "volume points | error                | " 
    << "prior / error        | lg2(prior/error) - 2 | runtime (s)"  
    << std::endl;

  std::cout << std::setprecision(14) << std::scientific;

  std::optional<long double> prior_error = std::nullopt;
  for (auto i = 0; i < 8; ++i) {

    std::vector<type::real_t> exact;
    poisson1d::analytical_blocked(exact, volume_points, blocks);

    std::vector<type::real_t> result;
    auto begin = timing::read();
    poisson1d::petsc_problem(result, volume_points, blocks);
    auto end = timing::read();
    auto secs = begin - end;

    long double error = 0;
    long double temp;
    for (std::size_t i = 0; i < result.size(); ++i) {
        temp = result[i] - exact[i];
        error += temp * temp;
    }
 
    error = std::sqrt(error  / (result.size() / 2)) ;

    std::cout << std::setw(13) 
      << exact.size() << " | " << error << " | ";
    if (prior_error) {
      std::cout << *prior_error / error << " | " 
        << std::log2(*prior_error / error) - 2. << " | " << secs;
    }
    std::cout << std::endl;

    prior_error = error;
    volume_points *= 2;
    volume_points += blocks - 1;
  }
}

void poisson_1d_hybrid_convergence_test(
  std::size_t volume_points, 
  std::size_t blocks) {

  /* Solving the Poisson Equation in 1D

    ⎧ u_xx π^2 sin(πx) = 0 , 0 < x < 1
    ⎨ u = 0    (Dirichlet) , x = 0 
    ⎩ u_x = -π (Neumann)   , x = 1 

  */

  std::cout << "[ Poisson 1D equation with arbitrary flux " 
    << "interfaces. ]" << std::endl;

  std::cout << "volume points | error    | δ error  | " 
  << "r error  | time to sol (s)" << std::endl;

  std::cout << std::setprecision(14) << std::scientific;

  std::optional<long double> prior_error = std::nullopt;
  for (auto i = 0; i < 8; ++i) {

    std::vector<type::real_t> exact;
    poisson1d::analytical_blocked(exact, volume_points, blocks);

    std::vector<type::real_t> result;
    auto begin = timing::read();
    poisson1d::petsc_hybridized_problem(result, volume_points, blocks);
    auto end = timing::read();
    auto secs = end - begin;

    long double error = 0;
    long double temp;
    for (std::size_t i = 0; i < result.size(); ++i) {
        temp = result[i] - exact[i];
        error += temp * temp;
    }
 
    error = std::sqrt(error  / (result.size() / 2)) ;

    std::cout << std::setprecision(2) << std::scientific;
    std::cout << std::setw(13) 
      << exact.size() << " | " << error << " | ";
    if (prior_error) {
      std::cout << *prior_error / error << " | " 
        << std::log2(*prior_error / error)  << " | ";
      std::cout << std::setprecision(4) << std::fixed;
      std::cout << secs;
    }
    std::cout << std::endl;

    prior_error = error;
    volume_points *= 2;
    volume_points += blocks - 1;
  }
}

void main_2d() {
  
  using namespace sbp_sat;

  nat_t bsz = 11;
  real_v u = {};

  block_v b = { 
    std::make_tuple(
      range_t(0,     1./3., bsz), 
      range_t(0,     1./3., bsz)),
    std::make_tuple(
      range_t(1./3., 2./3., bsz), 
      range_t(1./3., 2./3., bsz)),
    std::make_tuple(
      range_t(2./3., 1,     bsz), 
      range_t(2./3., 1.,    bsz))};

  domain_v d = {
    std::make_tuple(0., 1.), 
    std::make_tuple(0., 1.)};

  auto gw = [](long double x){return std::sin(x);};
  auto ge = [](long double x){return -std::sin(x);};
  auto gn = [](long double x){return -π * std::cos(π * x);};
  auto gs = [](long double x){return -π * std::cos(π * x);};

  boundary_vx2 g = {{
    std::make_tuple(range_t(0, 1./3., bsz) | gw, 1),
    std::make_tuple(range_t(0, 1./3., bsz) | ge, 1),
    std::make_tuple(range_t(0, 1./3., bsz) | gn, 2),
    std::make_tuple(range_t(0, 1./3., bsz) | gs, 2)
  }, {
    std::make_tuple(range_t(1./3., 2./3., bsz) | gw, 1),
    std::make_tuple(range_t(1./3., 2./3., bsz) | ge, 1),
    std::make_tuple(range_t(1./3., 2./3., bsz) | gn, 2),
    std::make_tuple(range_t(1./3., 2./3., bsz) | gs, 2)      
  }, {
    std::make_tuple(range_t(1./3., 2./3., bsz) | gw, 1),
    std::make_tuple(range_t(1./3., 2./3., bsz) | ge, 1),
    std::make_tuple(range_t(1./3., 2./3., bsz) | gn, 2),
    std::make_tuple(range_t(1./3., 2./3., bsz) | gs, 2)      
  }};

  boundary_tx2 gg = {
    std::make_tuple(range_t(0, 1., bsz) | gw, 1),
    std::make_tuple(range_t(0, 1., bsz) | ge, 1),
    std::make_tuple(range_t(0, 1., bsz) | gn, 2),
    std::make_tuple(range_t(0, 1., bsz) | gs, 2)};

  block_t bb = std::make_tuple(
    range_t(0, 1., bsz), 
    range_t(0, 1., bsz));

  // sbp_sat::x2::petsc_poisson(u, d, bb, gg);

  sbp_sat::x2::petsc_hybridized_poisson(u, d, b, g);

  
}
