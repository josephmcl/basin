/*
#include "timing.h"

int main(int argc, char **argv) {
    return 0;
}
*/

#include <hip/hip_runtime.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <numbers>
#include <chrono>


#include "rocalution.hpp"
#include "rocsolver.h"


#include "definitions.h"
#include "compute_boundary_solution.h"
#include "compute_sources.h"
#include "compute_g.h"
#include "compute_lambda_b.h"
// #include "compute_mg.h"
#include "compute_f_symbols.h"
#include "compute_rhs.h"
#include "cholesky.h"
#include "data.h"
#include "connect.h"

#include "omp.h"

#include <fstream>
#include <filesystem>

using real_t = type::real_t;

void write_chol(real_t *m, std::size_t sz, components &sbp, std::string key) {
  std::string dir = "../operator/";
  dir += "V_";
  dir += std::to_string(sbp.n * sbp.n * sbp.n_blocks_dim * sbp.n_blocks_dim); 
  dir += "_N_";
  dir += std::to_string(sbp.n); 
  dir += "_L_";
  dir += std::to_string(sbp.n_blocks_dim); 
  namespace fs = std::filesystem;
  fs::create_directories(dir);
  std::ofstream out;
  out.open(dir + "/" + key + ".cholesky_factors",
    std::ios::out | std::ios::binary);
  out.write(reinterpret_cast<const char*>(&m[0]), sizeof(real_t) * sz * sz);
  out.close();
  std::cout << dir + "/" + key + ".cholesky_factors " << std::endl;
}

bool load_factors(real_t *m, std::size_t sz, components &sbp, std::string key) {
  std::string dir = "../operator/";
  dir += "V_";
  dir += std::to_string(sbp.n * sbp.n * sbp.n_blocks_dim * sbp.n_blocks_dim); 
  dir += "_N_";
  dir += std::to_string(sbp.n); 
  dir += "_L_";
  dir += std::to_string(sbp.n_blocks_dim); 
  dir += "/";
  dir += key;
  dir +=".cholesky_factors";

  namespace fs = std::filesystem;
  const fs::path p{dir};
  if (fs::exists(p)) {
    auto ifs = std::ifstream(dir, std::ios::binary);
    ifs.read(
      reinterpret_cast<char*>(m), 
      sizeof(real_t) * sz * sz);
    ifs.close();
    double gb = (sizeof(real_t) * sz * sz) / 1e9;
    std::cout << "Loaded dense cholesky factors \"" << dir << "*\" (" 
    << std::fixed << gb << " GB)" << std::endl;
    return true;
  }
  else {
    return false;
  }
}

int main() {

  rocblas_handle rb_handle;
  rocblas_create_handle(&rb_handle);
  rocsparse_handle handle;
  rocsparse_create_handle(&handle);

  rocblas_status   rbstatus;
  rocsparse_status rsstatus;

  std::size_t vln = 125; //4; // 125 8 
  std::size_t eln = 8; //3;

  const std::size_t l_blocks = eln;
  const std::size_t n_blocks = l_blocks * l_blocks;
  const real_t span = 1. / static_cast<double>(l_blocks);
  const std::size_t n_points_x = vln;
  const std::size_t n = vln;

  const std::size_t block_size_dim = vln;
  const std::size_t block_count_dim = eln;

  vv<std::size_t> interfaces;
  std::size_t n_interfaces = make_connectivity(interfaces, l_blocks); 

  std::cout << "span " << span << std::endl;

  auto sbp = components{n, span};
  auto space = 0.5* (span/(n - 1));
  sbp.𝜏 = (2/ span) + (2 / (span * (space/span))); // * 10; 
    // std::cout << "TAU_VALUE " << sbp.TAU_VALUE << std::endl;
    // sbp.TAU_VALUE = (sbp.TAU_VALUE < 42.)? 42.: sbp.TAU_VALUE; // hard code these coeffs for now. 
  sbp.β = 1.;

  // int world_size;
  // MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  sbp.n_ranks = 1; // static_cast<std::size_t>(world_size);

  // int world_rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  sbp.rank = 0; // static_cast<std::size_t>(world_rank);

  sbp.n_blocks = n_blocks; // additional non sbp-sat info. but  
  sbp.n_interfaces = n_interfaces; // useful to have along. 
  sbp.n_blocks_dim = l_blocks;
  // sbp.n_threads = static_cast<std::size_t>(omp_get_max_threads());

  // alignment is the number of ranks that have to compute an extra
  // solution.
  std::size_t alignment = sbp.n_blocks % sbp.n_ranks;          
  sbp.rank_limit_u = (sbp.rank < alignment) 
      ? (sbp.n_blocks / sbp.n_ranks) + 1
      : (sbp.n_blocks / sbp.n_ranks);

  // slack adjusts the index the latter, smaller indices to align with
  // the larger, earlier indices. 
  std::size_t slack = (sbp.rank < alignment)? 0: std::min(sbp.rank, alignment);
  for (std::size_t i = 0; i != sbp.rank_limit_u; ++i) {
      sbp.rank_index_u.push_back(slack + sbp.rank_limit_u * sbp.rank + i);
  }

  auto gdata = data<real_t> {};

  // Generate ranges for the x and y of each block, given the number of 
  // blocks in each dimension.
  auto block_grid = range_t(0, 1., l_blocks + 1);
  std::vector<range_t> grids; 
  for (std::size_t i = 0; i != block_grid.size() - 1; ++i)  {
    grids.push_back(range_t(*block_grid[i], *block_grid[i + 1], n));
  }

  // Generate the solution at the boundary. This implementation generates 
    // vectors from range functions. 
    real_t *boundary_solution;
    compute_boundary_solution(  
        &boundary_solution, grids, {gdata.w, gdata.e, gdata.s, gdata.n}, {0., 1., 0., 1.});

    real_t *sources;
    compute_sources(&sources, grids, gdata.source_function);
  
  for (std::size_t i = 0; i < 144; ++i) {
  //  std::cout << sources[i] <<  " ";
  }
  std::cout << std::endl << std::endl;
    
  for (std::size_t i = 0; i < 48; ++i) {
  //  std::cout << boundary_solution[i] <<  " ";
  }
  std::cout << std::endl << std::endl;

  // compute_b(B, sbp);  

  csr_t M0, M1, M2, λA;

  auto B  = std::vector<csr_t>(8);
  auto F  = std::vector<csr_t>(4);
  auto FT = std::vector<csr_t>(4);

  load_operator(M0, "M0", block_size_dim, block_count_dim);
  load_operator(M1, "M1", block_size_dim, block_count_dim);
  load_operator(M2, "M2", block_size_dim, block_count_dim);
  load_operator(λA,  "T",  block_size_dim, block_count_dim);

  load_operator(B[0],  "B0",  block_size_dim, block_count_dim);
  load_operator(B[1],  "B1",  block_size_dim, block_count_dim);
  load_operator(B[2],  "B2",  block_size_dim, block_count_dim);
  load_operator(B[3],  "B3",  block_size_dim, block_count_dim);
  load_operator(B[4],  "B4",  block_size_dim, block_count_dim);
  load_operator(B[5],  "B5",  block_size_dim, block_count_dim);
  load_operator(B[6],  "B6",  block_size_dim, block_count_dim);
  load_operator(B[7],  "B7",  block_size_dim, block_count_dim);

  load_operator(F[0],  "F0",  block_size_dim, block_count_dim);
  load_operator(F[1],  "F1",  block_size_dim, block_count_dim);
  load_operator(F[2],  "F2",  block_size_dim, block_count_dim);
  load_operator(F[3],  "F3",  block_size_dim, block_count_dim);

  FT[0] = csr_t(F[0], true); /* Compute the transpose of each F   */
  FT[1] = csr_t(F[1], true); /* matrix.                           */
  FT[2] = csr_t(F[2], true);
  FT[3] = csr_t(F[3], true);

  vv<std::size_t> boundary_data_map;
  vv<std::size_t> boundary_order_map;
  make_boundary_maps(boundary_data_map, boundary_order_map, l_blocks);

  // Compute the hybrid system g terms. 
  real_t *g;
  compute_g(&g, B, boundary_solution, sources, boundary_order_map, 
    boundary_data_map, sbp);
  std::cout << "Computed g" << std::endl;

  for (std::size_t i = 0; i < sbp.n * sbp.n * sbp.rank_limit_u; ++i) {
  //  std::cout << g[i] <<  " ";
  }

  std::cout << std::endl;
  

  vv<std::size_t> F_symbols(n_blocks,           // rows 
  std::vector<std::size_t>(n_interfaces, 0)); // columns
  vv<std::size_t> FT_symbols(n_interfaces,      // rows
  std::vector<std::size_t>(n_blocks, 0));     // columns

  compute_f_symbols(F_symbols, FT_symbols, interfaces, sbp);
  std::cout << "Computed F symbolic matrix." << std::endl; 
  
  /*
  rocsparse_mat_descr descr;
  rocsparse_create_mat_descr(&descr);
  rocsparse_set_mat_type(descr, rocsparse_matrix_type_general);
  rocsparse_set_mat_index_base(descr, rocsparse_index_base_zero);
  
  rstatus = rocsparse_create_csr_descr(descr, M0.n, M0.m, M0.nnz(), 
    M0.row_index_data(), M0.col_index_data(), M0.val_data(), 
    rocsparse_indextype_i32, rocsparse_indextype_i32, 
    rocsparse_index_base_zero, rocsparse_datatype_f64_r);*/
 

  auto M = std::vector<double *>(3);
  double * λλ;
  M[0] = M0.dense();
  M[1] = M1.dense();
  M[2] = M2.dense();
  λλ = λA.dense();

  
  int info;

  if (!load_factors(M[0], M0.n, sbp, "M0")) {
    rbstatus = rocsolver_dpotf2(rb_handle, rocblas_fill_lower, M0.n, 
      M[0], M0.n, &info);
    hipDeviceSynchronize();

    if (rbstatus != 0) {
      std::cout << "rocsolver_dpotf2 (runtime error)" << std::endl;
    }
    else if (info != 0) {
      std::cout << "rocsolver_dpotf2 (factorization error)" << std::endl;
    }
    else {
      std::cout << "Factorized M[0]." << std::endl;
      write_chol(M[0], M0.n, sbp, "M0");
    }
  }

  if (!load_factors(M[1], M1.n, sbp, "M1")) {
    rbstatus = rocsolver_dpotf2(rb_handle, rocblas_fill_lower, M1.n, 
      M[1], M1.n, &info);
    hipDeviceSynchronize();

    if (rbstatus != 0) {
      std::cout << "rocsolver_dpotf2 (runtime error)" << std::endl;
    }
    else if (info != 0) {
      std::cout << "rocsolver_dpotf2 (factorization error)" << std::endl;
    }
    else {
      std::cout << "Factorized M[1]." << std::endl;
      write_chol(M[1], M1.n, sbp, "M1");
    }
  }
  
  if (!load_factors(M[2], M2.n, sbp, "M2")) {
    rbstatus = rocsolver_dpotf2(rb_handle, rocblas_fill_lower, M2.n, 
      M[2], M2.n, &info);
    hipDeviceSynchronize();  

    if (rbstatus != 0) {
      std::cout << "rocsolver_dpotf2 (runtime error)" << std::endl;
    }
    else if (info != 0) {
      std::cout << "rocsolver_dpotf2 (factorization error)" << std::endl;
    }
    else {
      std::cout << "Factorized M[2]." << std::endl;
      write_chol(M[2], M2.n, sbp, "M2");
    }
  }

  if (!load_factors(λλ, λA.n, sbp, "T")) {
    rbstatus = rocsolver_dpotf2(rb_handle, rocblas_fill_lower, λA.n, 
      λλ, λA.n, &info);
    hipDeviceSynchronize(); 
    
    if (rbstatus != 0) {
      std::cout << "rocsolver_dpotf2 (runtime error)" << std::endl;
    }
    else if (info != 0) {
      std::cout << "rocsolver_dpotf2 (factorization error)" << std::endl;
    }
    else {
      std::cout << "Factorized λA." << std::endl;
      write_chol(λλ, λA.n, sbp, "T");
    }
  }
  

  std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
  mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

  real_t *Mg, *gp, *mgp;
  std::size_t k, g_limit, mg_limit;
    
  g_limit = sbp.n * sbp.n * sbp.rank_limit_u;
  mg_limit = sbp.rank_limit_u;

  constexpr auto s = solver_t::dense_cholesky;

  sbp.rb_handle = rb_handle;
  
  //Mg = (real_t *) malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u);

  //#pragma omp for target teams distribute
  //for (std::size_t i = 0; i < g_limit; ++i) {
  // Mg[i] = g[i];
  //}

  //exit(-1);

  auto begin = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> dur;

  Mg = (real_t *) malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u);

  #pragma omp parallel for 
  for (std::size_t i = 0; i < g_limit; ++i) {
    Mg[i] = g[i];
  }

  begin = std::chrono::steady_clock::now();

  #pragma omp parallel for target teams distribute
  for (std::size_t i = 0; i < mg_limit; ++i) {

    /* Get the current thread index.                            */
    //auto td = omp_get_thread_num();

    /* Get index of the solution vector for the rank-local
        problem. Then, look up which M matrix is used to compute 
        the solution at this index.                              */
    auto ii = sbp.rank_index_u[i];
    k = mi[ii % sbp.n_blocks_dim];
    
    /* Look up the pointers to the correct                      */
    gp = &g[i * sbp.n * sbp.n];
    mgp = &Mg[i * sbp.n * sbp.n];
    
    /* Solve the local sub-problem.                             */
    if constexpr (s == solver_t::dense_cholesky) {
      auto error = cholesky::solve(M[k], mgp, sbp);
      if (error) {
          std::cout << "Choleksy solve failed code " 
              << *error << std::endl;
      }
    }
  }

  end = std::chrono::steady_clock::now();
  dur = end - begin;
  std::cout << "Solved M^-1 * x = g. " <<  dur << std::endl;

  // num_threads(sbp.n_threads)
  
  
  //for (std::size_t i = 0; i < g_limit; ++i) {
  //  std::cout << Mg[i] << std::endl;
  //}

  real_t *λb;
  std::size_t sz = sizeof(real_t) * sbp.n * sbp.n_interfaces;
  λb = (double *) malloc(sz);
  memset(λb, 0, sz);    

  // Compute λb
  begin = std::chrono::steady_clock::now();
  compute_lambda_b(λb, F, Mg, FT_symbols, sbp);
  end = std::chrono::steady_clock::now();
  dur = end - begin;
  std::cout << "Computed λb (gd - F^T * M \\ g). " << dur << std::endl;

  //for (std::size_t i = 0; i < sbp.n * sbp.n_interfaces; ++i) {
  //  std::cout << λb[i] <<  " ";
  //}
  //std::cout << std::endl;

  double *λ = nullptr;
  const std::size_t λ_limit = sbp.n_interfaces * sbp.n;
  λ = (double *) malloc(sizeof(double) * λ_limit);
  for (std::size_t i = 0; i < λ_limit; ++i) { λ[i] = λb[i]; }

  begin = std::chrono::steady_clock::now();
  if constexpr (s == solver_t::dense_cholesky) {
    
    auto error = cholesky::solve(λλ, λ, sbp, λA.n);
    if (error) {
        std::cout << "Choleksy solve failed code " 
            << *error << std::endl;
    }
    
  }
  end = std::chrono::steady_clock::now();
  dur = end - begin;
  std::cout << "Solved λA * λ = λb. " <<  dur << std::endl;

  //for (std::size_t i = 0; i < sbp.n_interfaces * sbp.n; ++i) {
  //  std::cout << λ[i] <<  " ";
  //}

  const std::size_t rhs_limit = sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u;
  real_t *rhs = (double *) malloc(rhs_limit);
  
  begin = std::chrono::steady_clock::now();
  compute_rhs(rhs, g, FT, λ, F_symbols, sbp);
  end = std::chrono::steady_clock::now();
  dur = end - begin;
  std::cout << "Computed rhs []. " <<  dur << std::endl;

  real_t *u; 
  // NOTE: rank_limit_max here so we can gather the same size on 
  //       everything.
  sz = sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u; 
  u = (double *) malloc(sz);
  memset(u, 0, sz);

  #pragma omp parallel default(none) shared(s, u, M, rhs, sbp, begin, end, mi, g_limit, mg_limit, M, dur) private(gp, mgp, k)
  { 
    /* Copy g into Mg if the solver requires it.                  */
    if constexpr (s == solver_t::dense_cholesky) {
      #pragma omp for 
      for (std::size_t i = 0; i < g_limit; ++i) {
        u[i] = rhs[i];
      }
    }

    #pragma omp for  
    for (std::size_t i = 0; i < g_limit; ++i) { u[i] = 0.; }

    /* */
    if constexpr (s == solver_t::dense_cholesky) {
      #pragma omp for 
      for (std::size_t i = 0; i < g_limit; ++i) { u[i] = rhs[i]; }
    }

      /* Retrieve timing at the beginning of Mx = g.                */
      #pragma omp single
      { 
        begin = std::chrono::steady_clock::now();
      }

      /* Compute Mx = g. For each sub-problem, grab the index of the 
         memory for M, x (i.e. mg), and g. Pass those indices into 
         the solver routine. */
      #pragma omp for 
      for (std::size_t i = 0; i < mg_limit; ++i) {

        /* Get the current thread index.                            */
        auto td = omp_get_thread_num();

        /* Get index of the solution vector for the rank-local
           problem. Then, look up which M matrix is used to compute 
           the solution at this index.                              */
        auto ii = sbp.rank_index_u[i];
        k = mi[ii % sbp.n_blocks_dim];
        
        /* Look up the pointers to the correct                      */
        gp = &rhs[i * sbp.n * sbp.n];
        mgp = &u[i * sbp.n * sbp.n];
        
        /* Solve the local sub-problem.                             */
        if constexpr (s == solver_t::dense_cholesky) {
          auto error = cholesky::solve_no_sync(M[k], mgp, sbp);
          if (error) {
              std::cout << "Choleksy solve failed code " 
                  << *error << std::endl;
          }
        }

      }

      /* Calculate timing for Mx = g.                               */
      #pragma omp single
      {
        // __itt_task_end(domain);
        hipDeviceSynchronize();
        end = std::chrono::steady_clock::now();
        //trace.push_back(end - begin);
        dur = end - begin;
        std::cout << "Computed u = M^-1 * rhs " <<  dur << std::endl;
      } 

      // #pragma omp single  
      // for (std::size_t i = 0; i < g_limit; ++i) {
      //   std::cout << Mg[i] << std::endl;
      // }

    } /* End OpenMP parallel region.                                */

}