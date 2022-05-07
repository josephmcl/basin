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
    auto to_real_t = [](nat_t n){return static_cast<real_t>(n);};

    void analytical(nat_t const size, real_t const left, 
                    real_t const right, vector_t &res) {
        auto h = (right - left) / static_cast<real_t>(size - 1);
        for (std::size_t i = 0; i < size; ++i) 
            res.push_back(std::sin(h * π * i));
        return;
    }

    domain_t const domain = std::make_tuple(0., 1.);
    boundary_t const boundary = std::make_tuple(0., 1, -π, 2);

    void problem(nat_t size, 
                 domain_t domain=domain, 
                 boundary_t boundary=boundary) {
        
        // Unpack tuples 
        auto [left, right] = domain;
        auto [left_data,  left_order, 
              right_data, right_order] = boundary;

        auto domain_range = range_t(left, right, size);   

        auto b = vector_t(size, 1.); 
        real_t i = 0;
        for (auto &e : b) {
            e = -π * π * std::sin(π * to_real_t(i++) / to_real_t(size - 1));
            
        }
        
        auto hi = numerical::operators::H_inverse(size, 2, left, right);

        b[size - 1] += right_data * hi[size - 1];

        for (auto &e : b) {
            std::cout << e << std::endl;
        }

        auto d2 = numerical::operators::d2(size, 2, left, right);
        
        d2.left_boundary(-40., 1);  // NOTE: -40 and 1 are sigma constants.
        d2.right_boundary(1., 2);
        //  d2.left_sat(1.);

        auto cg = solve::conj_grad {
            .size            = size,
            .matrix_action   = d2.product(),
            .inplace_vv_diff = linalg::ip_vv_diff<vector_t>,
            .vtv_prod        = linalg::vtv_prod<vector_t>,
            .vs_prod         = linalg::vs_prod<vector_t, real_t>,
            .inplace_vv_sum  = linalg::ip_vv_sum<vector_t>
        };

        auto x = vector_t(size, 0.); 
        type::real_t epsilon = 10e-10;
        cg(x, b, epsilon, size * size);

        // auto prod = d2.product();

        //prod(, x);

        (cg.converged)
            ? std::cout << "Converged." << std::endl
            : std::cout << "Failed to converge." << std::endl;
    }

  void petsc_problem(
    nat_t size, 
    nat_t local_problems=2,
    domain_t domain=domain, 
    boundary_t boundary=boundary) {
        
        long double β = 1.;
        long double σ1 = -40.;
        long double σ2 = 1;
        long double ϵ = 1.;

        // Unpack tuples 
        auto [left, right] = domain;
        auto [left_data,  left_order, 
              right_data, right_order] = boundary;

        auto domain_range = range_t(left, right, size);  

        auto local_domain_ranges = std::vector<range_t>();

        auto local_domain_size = (size + local_problems - 1) / local_problems;
        
        /* NOTE: The number of volume points differs from the number of
                 matrix points because the interface points are used in 
                 both problems. */

        auto volume_points = (local_domain_size - 1) *local_problems+1;
        auto matrix_points = local_domain_size * local_problems;

        if (volume_points != size) 
            std::cout << "NOTE: The problem was specified with " 
                << size << " points, but this does not evenly divide " 
                << "the number of local problems (" << local_problems
                <<  "), so it has been adjusted to " << volume_points
                << "." << std::endl;

        auto spacing = to_real_t(right - left) / (volume_points - 1);
        auto spacing_square = spacing * spacing;
        for (nat_t i = 0; i != local_problems; ++i) {
            
            auto begin = volume_points / local_problems * i * spacing;
            auto end   = volume_points / local_problems * (i + 1) * spacing;
            local_domain_ranges.push_back(
                range_t(begin, end, local_domain_size));
        }

    auto h = numerical::operators::H(local_domain_size, 2, 0, (right - left) / local_problems);
    auto hi = numerical::operators::H_inverse(local_domain_size, 2, 0, (right - left) / local_problems);
    vector_t bs = {(3./2.) / spacing, -2. / spacing, (1./2.) / spacing};

        
    /* Create a petsc object for the A matrix. */
    Mat A;
    MatCreateSeqAIJ(PETSC_COMM_SELF, matrix_points, matrix_points, 
        4, nullptr, &A);

    /* Set the second-order SBP D2 * H in A using D2 and H sized for 
       the local problem. */
    for (std::size_t i = 0; i < local_domain_ranges.size(); ++i) {
      auto local = local_domain_ranges[i];
      auto local_offset = i * local_domain_size;
      
      /* Initialize the first local skew row. */
      MatSetValue(A, local_offset, local_offset,  
         1. / spacing_square, ADD_VALUES);
      MatSetValue(A, local_offset, local_offset + 1, 
        -2. / spacing_square, ADD_VALUES);
      MatSetValue(A, local_offset, local_offset + 2, 
         1. / spacing_square, ADD_VALUES); 

      /* Initialize the final local skew row. */
      auto local_last_row = (i + 1) * local_domain_size;
      MatSetValue(A, local_last_row - 1, local_last_row - 3,  
         1. / spacing_square, ADD_VALUES);
      MatSetValue(A, local_last_row - 1, local_last_row - 2, 
        -2. / spacing_square, ADD_VALUES);
      MatSetValue(A, local_last_row - 1, local_last_row - 1, 
         1. / spacing_square, ADD_VALUES); 

      /* Initialize the interior local diagonal rows. */
      for (auto it = local.begin() + 1; it != local.end() - 1; ++it) {  
        auto global_row_index    = local_offset + it.index;
        auto global_column_index = local_offset + it.index - 1;

        MatSetValue(A, global_row_index, global_column_index,  
           1. / spacing_square, ADD_VALUES);
        MatSetValue(A, global_row_index, global_column_index + 1, 
          -2. / spacing_square, ADD_VALUES);
        MatSetValue(A, global_row_index, global_column_index + 2, 
           1. / spacing_square, ADD_VALUES);  
      }
    }

    /* Set the second-order SAT terms for the left hand first-order
       boundary: β * HI1 * BS' * e0 * e0' */
    MatSetValue(A, 0, 0, β * hi[0] * bs[0], ADD_VALUES);
    MatSetValue(A, 1, 0, β * hi[1] * bs[1], ADD_VALUES);
    MatSetValue(A, 2, 0, β * hi[2] * bs[2], ADD_VALUES);

    /* Set the boundary data for the left hand first-order boundary. */
    MatSetValue(A, 0, 0, σ1 * hi[0], ADD_VALUES);

    MatSetValue(A, matrix_points - 1, matrix_points - 3, 
      σ2 * hi[local_domain_size - 1] * bs[2], ADD_VALUES);
    MatSetValue(A, matrix_points - 1, matrix_points - 2, 
      σ2 * hi[local_domain_size - 1] * bs[1], ADD_VALUES);
    MatSetValue(A, matrix_points - 1, matrix_points - 1, 
      σ2 * hi[local_domain_size - 1] * bs[0], ADD_VALUES);


    for (std::size_t i = 0; i < local_domain_ranges.size() - 1; ++i) {
      auto local = local_domain_ranges[i];
      auto local_offset = (i + 1) * local_domain_size;
      auto n = local_domain_size;


      /* β * HI1 * BS' * en * en' */ 
      MatSetValue(A, local_offset - 3, local_offset - 1, 
        β * hi[n - 3] * bs[2], ADD_VALUES);
      MatSetValue(A, local_offset - 2, local_offset - 1, 
        β * hi[n - 2] * bs[1], ADD_VALUES);
      MatSetValue(A, local_offset - 1, local_offset - 1, 
        β * hi[n - 1] * bs[0], ADD_VALUES);

      /* ϵ * HI1 * BS  * en * en' */
      MatSetValue(A, local_offset - 1, local_offset - 3, 
        ϵ * hi[n - 1] * bs[2], ADD_VALUES);
      MatSetValue(A, local_offset - 1, local_offset - 2, 
        ϵ * hi[n - 1] * bs[1], ADD_VALUES);
      MatSetValue(A, local_offset - 1, local_offset - 1, 
        ϵ * hi[n - 1] * bs[0], ADD_VALUES);

      /* σ₁ * HI1 * en * en' */ 
      MatSetValue(A, local_offset - 1, local_offset - 1, 
        σ1 * hi[n - 1], ADD_VALUES);

      /* β * HI1 * BS' *e0 * e0' */
      MatSetValue(A, local_offset, local_offset, 
        β * hi[0] * bs[0], ADD_VALUES);
      MatSetValue(A, local_offset + 1, local_offset, 
        β * hi[1] * bs[1], ADD_VALUES);
      MatSetValue(A, local_offset + 2, local_offset, 
        β * hi[2] * bs[2], ADD_VALUES);

      /* ϵ * HI1 * BS * e0 * e0 */
      MatSetValue(A, local_offset, local_offset, 
        ϵ * hi[0] * bs[0], ADD_VALUES);
      MatSetValue(A, local_offset, local_offset + 1, 
        ϵ * hi[0] * bs[1], ADD_VALUES);
      MatSetValue(A, local_offset, local_offset + 2, 
        ϵ * hi[0] * bs[2], ADD_VALUES);

      /* σ₁ * HI1 * e0 * e0' */
      MatSetValue(A, local_offset, local_offset, 
        σ1 * hi[0], ADD_VALUES);

      /* -σ₁ * HI1 * en * e0' */
      MatSetValue(A, local_offset - 1, local_offset, 
        -σ1 * hi[n - 1], ADD_VALUES);

      /* -β * HI1 * BS' * en * e0' */
      MatSetValue(A, local_offset - 3, local_offset, 
        -β * hi[n - 3] * bs[2], ADD_VALUES);
      MatSetValue(A, local_offset - 2, local_offset, 
        -β * hi[n - 2] * bs[1], ADD_VALUES);
      MatSetValue(A, local_offset - 1, local_offset, 
        -β * hi[n - 1] * bs[0], ADD_VALUES);


      /* -β * HI1 * BS' * e0 * en' */
      MatSetValue(A, local_offset + 2, local_offset - 1, 
        -β * hi[2] * bs[2], ADD_VALUES);
      MatSetValue(A, local_offset + 1, local_offset - 1, 
        -β * hi[1] * bs[1], ADD_VALUES);
      MatSetValue(A, local_offset, local_offset, - 1 
        -β * hi[0] * bs[0], ADD_VALUES);


      /* -σ₁ * HI1 * e0 * en' */
      MatSetValue(A, local_offset, local_offset - 1, 
        -σ1 * hi[0], ADD_VALUES);

      /* ϵ * HI1 * en * e0' * BS */
       MatSetValue(A, local_offset - 1, local_offset, 
        ϵ * hi[n - 1] * bs[0], ADD_VALUES);
      MatSetValue(A, local_offset - 1, local_offset + 1, 
        ϵ * hi[n - 1] * bs[1], ADD_VALUES);
      MatSetValue(A, local_offset - 1, local_offset + 2, 
        ϵ * hi[n - 1] * bs[2], ADD_VALUES);

      /* ϵ * HI1 * e0 * en' * BS */
       MatSetValue(A, local_offset, local_offset - 1, 
        ϵ * hi[0] * bs[0], ADD_VALUES);
      MatSetValue(A, local_offset, local_offset - 2, 
        ϵ * hi[0] * bs[1], ADD_VALUES);
      MatSetValue(A, local_offset, local_offset - 3, 
        ϵ * hi[0] * bs[2], ADD_VALUES);

      
    }
    
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    MatView(A,PETSC_VIEWER_STDOUT_WORLD);


    Vec b;
    VecCreateSeq(PETSC_COMM_SELF, matrix_points, &b);

    for (std::size_t i = 0; i < local_domain_ranges.size(); ++i) {
      auto local = local_domain_ranges[i];
      auto local_offset = i * local_domain_size;
      for (auto it = local.begin(); it != local.end(); ++it) {  
        auto index = local_offset + it.index;
        VecSetValue(b, index, -π * π * std::sin(*it * π), ADD_VALUES);
      }
    }

    Vec x;
    VecCreateSeq(PETSC_COMM_SELF, matrix_points, &x);

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetUp(ksp);
    KSPSolve(ksp, b, x);
    VecView(x, PETSC_VIEWER_STDOUT_WORLD);

    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&x);

  }
};
