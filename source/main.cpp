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

#include "petscsys.h"

// #include "timing.h"

#include "hybrid_sbp_sat_2d.h"

void poisson_1d_convergence_test(
  std::size_t volume_points, 
  std::size_t blocks);

void poisson_1d_hybrid_convergence_test(
  std::size_t volume_points, 
  std::size_t blocks);

void main_2d();

int main (int argc, char *argv[]) {

  // Initialize petsc, among other libraries.
  if (!infrastructure::initialize(argc, argv)) exit(-1);; { 

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
    // auto begin = timing::read();
    poisson1d::petsc_problem(result, volume_points, blocks);
    // auto end = timing::read();
    // auto secs = begin - end;

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
      //std::cout << *prior_error / error << " | " 
      //  << std::log2(*prior_error / error) - 2. << " | " << secs;
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
    //auto begin = timing::read();
    poisson1d::petsc_hybridized_problem(result, volume_points, blocks);
    //auto end = timing::read();
    // auto secs = end - begin;

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
      // std::cout << secs;
    }
    std::cout << std::endl;

    prior_error = error;
    volume_points *= 2;
    volume_points += blocks - 1;
  }
}

void main_2d() {
  
  // 83160 7560 840 280, 3

  sbp_sat::x2::petsc_hybridized_poisson(50, 3);
  
}
