#include "poisson_2d.h"
#include <ctime>
#include <math.h>
#include <iomanip>
#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>

//#include "ittnotify.h"

#include "mkl_spblas.h"

void print_csr(sparse_matrix_t *m, std::size_t sz, components &sbp, std::string key) {

  sparse_matrix_t mm;
  auto s = mkl_sparse_convert_csr(*m, SPARSE_OPERATION_NON_TRANSPOSE, &mm);
  mkl_sparse_status(s);

  int ssz = static_cast<int>(sz);

  std::string dir = "./operator/";
  dir += "V_";
  dir += std::to_string(sbp.n * sbp.n * sbp.n_blocks_dim * sbp.n_blocks_dim); 
  dir += "_N_";
  dir += std::to_string(sbp.n); 
  dir += "_L_";
  dir += std::to_string(sbp.n_blocks_dim); 

  namespace fs = std::filesystem;
  
  fs::create_directories(dir);

  int *rs, *re, *ci; 
  double *v;
  sparse_index_base_t id = SPARSE_INDEX_BASE_ZERO;
  s = mkl_sparse_d_export_csr(mm, &id, &ssz, &ssz, &rs, &re, &ci, &v);
  mkl_sparse_status(s);

  //out.write( reinterpret_cast<const char*>( &f ), sizeof( float ));
  
  std::ofstream out;
  out.open( dir + "/" + key + ".row", std::ios::out | std::ios::binary);
  out.write( reinterpret_cast<const char*>(&rs[0]), sizeof(int) * (sz + 1));
  out.close();
  //std::ifstream ifs( dir + "/" + key + ".row", std::ios::binary);
  //int read[100];
  //ifs.read( reinterpret_cast<char*>(&read), sizeof(int) * (sz + 1));
  //for (std::size_t i = 0; i < sz + 1; ++i) {
  //  std::cout << read[i] << ' ';
  //}
  out.open(dir + "/" + key + ".col", std::ios::out | std::ios::binary);
  out.write(reinterpret_cast<const char*>(&ci[0]), sizeof(int) * (rs[sz]));
  out.close();
  out.open(dir + "/" + key + ".val", std::ios::out | std::ios::binary);
  out.write(reinterpret_cast<const char*>(&v[0]), sizeof(double) * (rs[sz]));
  out.close();

  out.open(dir + "/" + key + ".rsize", std::ios::out | std::ios::binary);
  out.write(reinterpret_cast<const char*>(&sz), sizeof(int));
  out.close();

  out.open(dir + "/" + key + ".csize", std::ios::out | std::ios::binary);
  out.write(reinterpret_cast<const char*>(&rs[sz]), sizeof(int));
  out.close();


  //for (std::size_t i = 0; i < rs[sz]; ++i) {
  //  std::cout << std::setprecision(18) << v[i] << " ";
  //}
  //std::cout << std::endl << std::endl; 

  //std::ifstream ifs(dir + "/" + key + ".val", std::ios::binary);
  //double read[10000];
  //ifs.read( reinterpret_cast<char*>(&read), sizeof(double) * rs[sz]);
  //for (std::size_t i = 0; i < rs[sz]; ++i) {
  //  std::cout << read[i] << ' ';
  //}
  //std::cout << std::endl << std::endl; 
}

void poisson_2d::problem(std::size_t vln, std::size_t eln) {

    timing::init();

    std::cout << "OMP max threads: " <<  omp_get_max_threads() << std::endl;

    // omp_set_nested(1);

    const std::size_t l_blocks = eln;
    const std::size_t n_blocks = l_blocks * l_blocks;
    auto span = 1. / static_cast<double>(l_blocks);
    

    auto n_points_x = vln;
    auto n = vln;

    std::cout << "span " << span << std::endl;


    auto sbp = components{n, span};

    auto space = 0.5* (span/(n - 1));
    sbp.TAU_VALUE = (2/ span) + (2 / (span * (space/span))); // * 10; 
    // std::cout << "TAU_VALUE " << sbp.TAU_VALUE << std::endl;
    // sbp.TAU_VALUE = (sbp.TAU_VALUE < 42.)? 42.: sbp.TAU_VALUE; // hard code these coeffs for now. 
    sbp.BETA_VALUE = 1.;
  
    std::cout << "TAU = " << sbp.TAU_VALUE << std::endl;

  
    auto gw = [](real_t x, real_t y){return std::sin(PI_VALUE * x + PI_VALUE * y);};
    auto ge = [](real_t x, real_t y){return std::sin(PI_VALUE * x + PI_VALUE * y);};
    auto gs = [](real_t x, real_t y){return -PI_VALUE * std::cos(PI_VALUE * x + PI_VALUE * y);};
    auto gn = [](real_t x, real_t y){return PI_VALUE * std::cos(PI_VALUE * x + PI_VALUE * y);};
    
    auto polar_to_cartesian = [](real_t r, real_t theta) {
      return std::make_tuple(r * cos(theta), r * sin(theta));
    };


    // Multiply by 2 here because u_xx == u_yy
    
    auto source_function = [](real_t x, real_t y) {
        return -2. * PI_VALUE * PI_VALUE * sin(PI_VALUE * x + PI_VALUE * y);};

    vv<std::size_t> interfaces;
    std::size_t n_interfaces = make_connectivity(interfaces, l_blocks); 

    sbp.n_blocks = n_blocks; // additional non sbp-sat info. but  
    sbp.n_interfaces = n_interfaces; // useful to have along. 
    sbp.n_blocks_dim = l_blocks;
    sbp.n_threads = static_cast<std::size_t>(omp_get_max_threads());
    

    //
    // Compute MPI parameters.
    //

    

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    sbp.n_ranks = static_cast<std::size_t>(world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    sbp.rank = static_cast<std::size_t>(world_rank);

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

    //
    // End MPI parameters.
    //

    // mkl_set_num_threads(sbp.n_threads);

    if (sbp.n_blocks == 1) {
      std::cout << "SINGLE BLOCK CASE" << std::endl 
        << " | Problem size: " << n_points_x << " x " << n_points_x 
        << std::endl;
      single(sbp);
      return;
    }

    // Compute flux components used by b. 
    std::vector<sparse_matrix_t> B;
    compute_b(B, sbp);  

    // Generate ranges for the x and y of each block, given the number of 
    // blocks in each dimension.
    auto block_grid = range_t(0, 1., sbp.n_blocks_dim + 1);
    std::vector<range_t> grids; 
    for (std::size_t i = 0; i != block_grid.size() - 1; ++i)  {
        real_t from = *block_grid[i];
        real_t to = *block_grid[i + 1];
        std::cout << from << " - " << to << std::endl;
        grids.push_back(range_t(from, to, n));
    }

    // Generate the solution at the boundary. This implementation generates 
    // vectors from range functions. 
    real_t *boundary_solution;
    compute_boundary_solution(  
        &boundary_solution, grids, {gw, ge, gs, gn}, {0., 1., 0., 1.});

    real_t *sources;
    compute_sources(&sources, grids, source_function);

    /*
    for (std::size_t i = 0; i < sbp.n_blocks_dim; ++i) {
      for (std::size_t j = 0; j < sbp.n; ++j) {
        for (std::size_t k = 0; k < sbp.n_blocks_dim; ++k) {
          for (std::size_t l = 0; l < sbp.n; ++l) {
            auto index = 
                (i * sbp.n * sbp.n * sbp.n_blocks_dim) 
              + (k * sbp.n * sbp.n) 
              + (j * sbp.n) 
              + l;
            std::cout << sources[index] << " ";
          }
        }
        std::cout << std::endl;
      }
    }
    */
    
    vv<std::size_t> boundary_data_map;
    vv<std::size_t> boundary_order_map;
    make_boundary_maps(boundary_data_map, boundary_order_map, l_blocks);


    // Compute the hybrid system g terms. 
    real_t *g;
    compute_g(&g, B, boundary_solution, sources, boundary_order_map, 
    boundary_data_map, sbp); 

    vv<std::size_t> F_symbols(n_blocks,           // rows 
    std::vector<std::size_t>(n_interfaces, 0)); // columns
    vv<std::size_t> FT_symbols(n_interfaces,      // rows
    std::vector<std::size_t>(n_blocks, 0));     // columns

    compute_f_symbols(F_symbols, FT_symbols, interfaces, sbp);

    if (sbp.rank == 0) {
      std::cout << "Square hybrid specs:" << std::endl 
      << " | local problem size: " << n_points_x << " x " << n_points_x 
      << std::endl << " | span: " << span 
      << std::endl << " | total blocks: " << l_blocks << " x " << l_blocks 
        << " (" << n_blocks << ") "
      << std::endl << " | total interfaces: " << n_interfaces 
      << std::endl << " | threads: " << sbp.n_threads 
      << std::endl << " | ranks: " << sbp.n_ranks
      << std::endl;
    }

    auto trace = std::vector<double>();
    
    auto M = vv<sparse_matrix_t>();
    M.resize(sbp.n_threads + 2);
    for (auto &e: M) {
      e.resize(3);
    }
    auto begin = timing::read();
    make_m(&M[0][0], sbp, {1, 1, 2, 1});
    make_m(&M[0][1], sbp, {1, 1, 1, 1});
    make_m(&M[0][2], sbp, {1, 1, 1, 2});
    auto end = timing::read();
    trace.push_back(end - begin);

    print_csr(&M[0][0], n_points_x * n_points_x, sbp, "M0");
    print_csr(&M[0][1], n_points_x * n_points_x, sbp, "M1");
    print_csr(&M[0][2], n_points_x * n_points_x, sbp, "M2");

    // logging::out << std::setw(14) << std::fixed << end - begin
    //  << " s # " << "Computed M." << std::endl;
    sparse_status_t status;
    matrix_descr dc;
    dc.type = SPARSE_MATRIX_TYPE_GENERAL;
    for (std::size_t i = 1; i != M.size(); ++i) {
      status = mkl_sparse_copy(M[0][0], dc, &M[i][0]);
      mkl_sparse_status(status);
      status = mkl_sparse_copy(M[0][1], dc, &M[i][1]);
      mkl_sparse_status(status);
      status = mkl_sparse_copy(M[0][2], dc, &M[i][2]);
      mkl_sparse_status(status);
    }

//#ifdef BASIN_CHOLESKY
//    auto constexpr s = solver::dense_cholesky;
//#else
    auto constexpr s = solver::sparse_qr;
//#endif

    std::cout << "Active solver - ";
    if constexpr (s == solver::dense_cholesky) {
      std::cout << "LAPACK dense Cholesky" << std::endl;
    }
    else if constexpr (s == solver::rfp_cholesky) {
      std::cout << "LAPACK RFP Cholesky" << std::endl;
    }
    else if constexpr (s == solver::sparse_qr) {
      std::cout << "Intel MKL sparse QR" << std::endl;
    }
    else if constexpr (s == solver::pardiso) {
      std::cout << "PARDISO" << std::endl;
    }

    std::vector<dense_matrix_z> 
      M_dense_cholesky, M_rfp_cholesky; 
    M_dense_cholesky.resize(3);
    M_rfp_cholesky.resize(3);

    if constexpr (s == solver::dense_cholesky) {
      for (std::size_t i = 0; i != 3; ++i) 
        M_dense_cholesky[i] = factor<solver::dense_cholesky>::ize(
          sbp, M[sbp.n_threads][i]);
    }
    else if constexpr (s == solver::rfp_cholesky) {
      for (std::size_t i = 0; i != 3; ++i) 
        M_rfp_cholesky[i] = factor<solver::dense_cholesky>::ize(
          sbp, M[sbp.n_threads][i]);
    }
    else if constexpr (s == solver::sparse_qr) {
      for (std::size_t i = 0; i != M.size(); ++i) {
        for (std::size_t j = 0; j != M[i].size(); ++j) {
          factor<solver::sparse_qr>::ize(sbp, M[i][j]);
        }
      }
    }
    

    // Compute F components.
    auto Fsparse = std::vector<sparse_matrix_t>(4);
    auto Fdense = std::vector<real_t *>(4);
    begin = timing::read();
    compute_f(Fsparse, Fdense, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin 
    // << " s # " << "Computed F." << std::endl;

    std::vector<real_t *> MF;
    MF.resize(M[0].size() * Fdense.size());
    for (std::size_t index = 0; index != MF.size(); ++index) {
        MF[index] = (real_t *) mkl_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.n, 64);
        memset(MF[index], 0, sizeof(real_t) * sbp.n * sbp.n * sbp.n);
    }
    
    if constexpr(s != solver::sparse_qr) {
      for (std::size_t i = 0; i < M[0].size(); ++i) {
        for (std::size_t j = 0; j < Fdense.size(); ++j) {
          for (std::size_t k = 0; k < sbp.n * sbp.n * sbp.n; ++k) {
            MF[i * Fdense.size() + j][k] = Fdense[j][k];
          }
        }
      }
    }

    // Compute solve of MX = F.
    begin = timing::read();
    if constexpr (s == solver::sparse_qr) {
      mf<solver::sparse_qr>::solve(sbp, MF, M, Fdense);
    }
    else if constexpr (s == solver::dense_cholesky) {
      mf<solver::dense_cholesky>::solve(sbp, MF, M_dense_cholesky);
    }
    else if constexpr (s == solver::rfp_cholesky) {
      mf<solver::rfp_cholesky>::solve(sbp, MF, M_rfp_cholesky);
    }
    end = timing::read();
    trace.push_back(end - begin);
    logging::out << "Computed MX=F." << std::endl;

    /*
    for (std::size_t k = 0; k < sbp.n * sbp.n * sbp.n; ++k) {
      std::cout << MF[0][k]  << " ";
    }
    std::cout << std::endl;
    */

    // Setup D matrix.
    sparse_matrix_t D;
    compute_d(&D, sbp, interfaces);
    // std::cout << "computed d " << std::endl;

    // Setup interface list.
    vv<std::size_t> lambda_indices;
    make_interface_list(lambda_indices, F_symbols, FT_symbols, sbp);
        
    double *Lambda_A = nullptr;
    MKL_INT *piv = nullptr;
    if (sbp.rank == 0) {
      // Allocate memory for Lambda_A.
      Lambda_A = (double *) mkl_malloc(
          sizeof(double) * sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces, 64);
      for (std::size_t i = 0; i != sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces; ++i) {
        Lambda_A[i] = 0;
      }
      
      // Allocate pivot table memory.
      //piv = (MKL_INT *) mkl_malloc(
      //  sizeof(MKL_INT) * sbp.n * sbp.n_interfaces, 64);
      //memset(piv, 0, sizeof(MKL_INT) * sbp.n * sbp.n_interfaces);
      
    }

    sparse_matrix_t Lambda_A_sparse;
    if (sbp.rank == 0) {
      
      // Compute Lambda_A.
      begin = timing::read();
      compute_lambda_a(Lambda_A, &D, Fsparse, MF, F_symbols, FT_symbols, sbp);
      end = timing::read();
      trace.push_back(end - begin);

      const int lsize = sbp.n_interfaces * sbp.n;

      double *av;
      int nnz, *ai, *aj;
      if constexpr (s == solver::sparse_qr) {
        nnz = ((sbp.n * sbp.n)) * (
          (16 * sbp.n_blocks) - (28 * sbp.n_blocks_dim) + 8);
        MKL_INT job[6] = {0, 0, 0, 2, nnz, 1};      
        av = (double *) mkl_malloc(sizeof(double) * nnz, 64);
        aj = (int *) mkl_malloc(sizeof(int) * nnz, 64);
        ai = (int *) mkl_malloc(sizeof(int) * lsize + 1, 64);
        int info;
        mkl_ddnscsr((int *)&job , &lsize, &lsize, *&Lambda_A, &lsize, 
          av, aj, ai, &info);
        
        status = mkl_sparse_d_create_csr(
              &Lambda_A_sparse, 
              SPARSE_INDEX_BASE_ZERO,
              lsize, lsize,
              &ai[0],
              &ai[1],
              &aj[0],
              &av[0]);
        mkl_sparse_status(status);
      }

      print_csr(&Lambda_A_sparse, lsize, sbp, "T");

            
      /* Factorize global lambda matrix.                            * 
       *                                                            */
      begin = timing::read();
      if constexpr (s == solver::dense_cholesky) {
        auto error = cholesky::factor(
          Lambda_A, sbp, sbp.n_interfaces * sbp.n);
        if (error) {
            std::cout << "Global cholesky factor failed "
                << "with code " << *error << std::endl; 
        } 
      }
      else if constexpr (s == solver::rfp_cholesky) {
        // Lambda_A = factor<solver::rfp_cholesky>::ize(
        //   sbp, Lambda_A);
      }
      else if constexpr (s == solver::sparse_qr) {
        factor<solver::sparse_qr>::ize(sbp, Lambda_A_sparse);
      }
      end = timing::read();

      logging::out << "Global matrix factorized." << std::endl;
      trace.push_back(end - begin);
    }
    else {
      trace.push_back(0);
      trace.push_back(0);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    /* Compute Mx = g directly via M factors.                       */ 
                                                                      
    
    /* mi takes the layout [0, 1, ..., 1, 2] and is of length l, i.e.,
       n_blocks_dim, or the number of sub-problems in one dimension of
       the problem. This is used to find the right sub-problem matrix
       from an artibtrary sub-problem index.                        */
    std::vector<std::size_t> mi(sbp.n_blocks_dim, 1);
    mi[0] = 0; mi[sbp.n_blocks_dim - 1] = 2;

    real_t *Mg, *gp, *mgp;
    std::size_t k, g_limit, mg_limit;
    
    g_limit = sbp.n * sbp.n * sbp.rank_limit_u;
    mg_limit = sbp.rank_limit_u;

    //__itt_domain* domain;
    #pragma omp parallel num_threads(sbp.n_threads) default(none) shared(s, Mg, M, g, sbp, begin, end, trace, status, mi, g_limit, mg_limit, M_dense_cholesky) private(gp, mgp, k)
    { 
      #pragma omp single
      {
        Mg = mg<s>::allocate(sbp);
      }

      /* Copy g into Mg if the solver requires it.                  */
      if constexpr (s == solver::dense_cholesky) {
        #pragma omp for 
        for (std::size_t i = 0; i < g_limit; ++i) {
          Mg[i] = g[i];
        }
      }

      /* Call each solver once.                                     */
      #pragma omp for 
      for (std::size_t i = 0; i < mg_limit; ++i) {
        auto td = omp_get_thread_num();
        auto ii = sbp.rank_index_u[i];
        gp = &g[i * sbp.n * sbp.n];
        mgp = &Mg[i * sbp.n * sbp.n];
        // rank_index_u is the global index, thats what we modulo. 
        k = mi[ii % sbp.n_blocks_dim];
        
        if constexpr (s == solver::dense_cholesky) {
          auto error = cholesky::solve(M_dense_cholesky[k], mgp, sbp);
          if (error) {
              std::cout << "M^(-1) g choleksy solve failed code " 
                  << *error << std::endl;
          }
        }
        else if constexpr (s == solver::sparse_qr) {
          status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
          mkl_sparse_status(status);
        }
      }

      #pragma omp for  
      for (std::size_t i = 0; i < g_limit; ++i) { Mg[i] = 0.; }

      /* */
      if constexpr (s == solver::dense_cholesky) {
        #pragma omp for 
        for (std::size_t i = 0; i < g_limit; ++i) { Mg[i] = g[i]; }
      }

      /* Retrieve timing at the beginning of Mx = g.                */
      #pragma omp single
      { 
        begin = timing::read();
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
        gp = &g[i * sbp.n * sbp.n];
        mgp = &Mg[i * sbp.n * sbp.n];
        
        /* Solve the local sub-problem.                             */
        if constexpr (s == solver::dense_cholesky) {
          auto error = cholesky::solve(M_dense_cholesky[k], mgp, sbp);
          if (error) {
              std::cout << "Choleksy solve failed code " 
                  << *error << std::endl;
          }
        }
        else if constexpr (s == solver::sparse_qr) {
          status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
          mkl_sparse_status(status);
        }

      }

      /* Calculate timing for Mx = g.                               */
      #pragma omp single
      {
        // __itt_task_end(domain);
        end = timing::read();
        trace.push_back(end - begin);
        logging::out << "Computed x = M^-1 * g." << std::endl;
      } 

      // #pragma omp single  
      // for (std::size_t i = 0; i < g_limit; ++i) {
      //   std::cout << Mg[i] << std::endl;
      // }

    } /* End OpenMP parallel region.                                */

    std::size_t sz = sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u;

    

    real_t *LAMBDAb;
    sz = sizeof(real_t) * sbp.n * sbp.n_interfaces;
    LAMBDAb = (double *) mkl_malloc(sz, 64);
    memset(LAMBDAb, 0, sz);    

    // Compute LAMBDAb
    begin = timing::read();
    compute_lambda_b(LAMBDAb, Fsparse, Mg, FT_symbols, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin 
    // << " s # " << 
    logging::out << "Computed LAMBDAb (gd - FT * M \\ g)." << std::endl;

    real_t *LAMBDAb_reduced = nullptr;
    sz = sizeof(real_t) * sbp.n * sbp.n_interfaces;
    LAMBDAb_reduced = (double *) mkl_malloc(sz, 64);

    int err = 0;
    sz = sbp.n * sbp.n_interfaces;

    begin = timing::read();
    err = MPI_Reduce(LAMBDAb, LAMBDAb_reduced, sz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    end = timing::read();
    
    trace.push_back(end - begin);

    if (err != 0 ) {
      std::cout << sbp.rank << " MPI REDUCE SUM ERROR CODE: " << 
        err << std::endl;
    }


    double *lamu = nullptr;
    lamu = (double *) mkl_malloc(sizeof(double) * sbp.n_interfaces * sbp.n, 64);
    
    /*
    double *ones = nullptr;
    ones = (double *) mkl_malloc(sizeof(double) * sbp.n_interfaces * sbp.n, 64);
    for (std::size_t i = 0; i < sbp.n_interfaces * sbp.n; ++i) {
        ones[i] = 1.;
      }
    */
    
    if (sbp.rank == 0) {
      
      /*
      begin = timing::read();
      status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, Lambda_A_sparse, nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, lamu , sbp.n * sbp.n_interfaces, 
            ones, sbp.n * sbp.n_interfaces);
      end = timing::read();

      std::cout << "lambda time 1: " << end - begin << std::endl;
      mkl_sparse_status(status);
      for (std::size_t i = 0; i < sbp.n_interfaces * sbp.n; ++i) {
        lamu[i] = 0.;
      }
      */

      if constexpr (s == solver::dense_cholesky) {
        for (std::size_t i = 0; i < sbp.n * sbp.n_interfaces; ++i) {
          lamu[i] = LAMBDAb[i];
        }  
      }

      begin = timing::read();
      if constexpr (s == solver::sparse_qr) {
        status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, Lambda_A_sparse, nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, lamu , sbp.n * sbp.n_interfaces, 
            LAMBDAb, sbp.n * sbp.n_interfaces);
        mkl_sparse_status(status);
      }
      else if constexpr (s == solver::dense_cholesky) {
        auto error = cholesky::solve(Lambda_A, lamu, sbp, 
          sbp.n * sbp.n_interfaces);
        if (error) {
          std::cout << "Global choleksy solve failed code " 
            << *error << std::endl;
        }
      }
      end = timing::read();

      logging::out << "Lambda solve completed." << std::endl;
      trace.push_back(end - begin);
    }
    else {
      trace.push_back(0);
    }

    //for (int i = 0; i < 20; ++i) {
    //  std::cout << altlu[i] - LAMBDAb[i] << " ";
    //}
    //std::cout << std::endl; 

    sz = sizeof(real_t) * sbp.n * sbp.n_interfaces;
    if (sbp.rank == 0) {
      std::memcpy(LAMBDAb, LAMBDAb_reduced, sz);
      mkl_free(LAMBDAb_reduced);
    }
    else {
      mkl_free(LAMBDAb_reduced);
    }

    begin = timing::read();
    sz = sbp.n * sbp.n_interfaces;
    MPI_Bcast(LAMBDAb, sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    end = timing::read();

    trace.push_back(end - begin);

    // Allocated right-hand side memory. 
    real_t *rhs;
    sz = sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u;
    rhs = (double *) mkl_malloc(sz, 64);
    memset(rhs, 0, sz);

    // Compute right-hand side. b = g - F * LAMBDA.
    begin = timing::read();
    compute_rhs(rhs, g, Fsparse, lamu, F_symbols, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin 
    //   << " s # 
    logging::out <<  "Computed b (g - F * LAMBDA)." << std::endl;

    //for (std::size_t i = 0; i < sbp.n_blocks * sbp.n * sbp.n; ++i) {
      //if (isnan(rhs[i])) {
    //    std::cout << rhs[i] << " " << std::endl;
      //}
    //}

    // Free g, F (sparse), and LAMBDA.
    mkl_free(g);
    for (auto &e: Fsparse) 
      mkl_sparse_destroy(e);
    mkl_free(LAMBDAb);
    mkl_free(lamu);

    // Allocate solution memory. 
    real_t *u; 
    // NOTE: rank_limit_max here so we can gather the same size on 
    //       everything.
    sz = sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u; 
    u = (double *) mkl_malloc(sz, 64);
    memset(u, 0, sz);
    
    if constexpr(s != solver::sparse_qr) {
      for (std::size_t i = 0; i < sbp.n * sbp.n * sbp.rank_limit_u; ++i) {
        u[i] = rhs[i];
      }
    }

    // std::cout << "banana" << std::endl;

    // Compute solution. u = M \ b.
    /*
    begin = timing::read();
    if constexpr (s == solver::sparse_qr) {
      mrhs<solver::sparse_qr>::solve(sbp, u, M, rhs);
    }
    else if constexpr (s == solver::dense_cholesky) {
      mrhs<solver::dense_cholesky>::solve(sbp, u, M_dense_cholesky);
    }
    else if constexpr (s == solver::rfp_cholesky) {
      mrhs<solver::rfp_cholesky>::solve(sbp, u, M_rfp_cholesky);
    }
    end = timing::read();

    trace.push_back(end - begin);
    */

    #pragma omp parallel num_threads(sbp.n_threads) default(none) shared(s, u, M, rhs, sbp, begin, end, trace, status, mi, g_limit, mg_limit, M_dense_cholesky) private(gp, mgp, k)
    { 
      /* Copy g into Mg if the solver requires it.                  */
      if constexpr (s == solver::dense_cholesky) {
        #pragma omp for 
        for (std::size_t i = 0; i < g_limit; ++i) {
          u[i] = rhs[i];
        }
      }

      /* Call each solver once.                                     */
      #pragma omp for 
      for (std::size_t i = 0; i < mg_limit; ++i) {
        auto td = omp_get_thread_num();
        auto ii = sbp.rank_index_u[i];
        gp = &rhs[i * sbp.n * sbp.n];
        mgp = &u[i * sbp.n * sbp.n];
        // rank_index_u is the global index, thats what we modulo. 
        k = mi[ii % sbp.n_blocks_dim];
        
        if constexpr (s == solver::dense_cholesky) {
          auto error = cholesky::solve(M_dense_cholesky[k], mgp, sbp);
          if (error) {
              std::cout << "M^(-1) rhs choleksy solve failed code " 
                  << *error << std::endl;
          }
        }
        else if constexpr (s == solver::sparse_qr) {
          status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
          mkl_sparse_status(status);
        }
      }

      #pragma omp for  
      for (std::size_t i = 0; i < g_limit; ++i) { u[i] = 0.; }

      /* */
      if constexpr (s == solver::dense_cholesky) {
        #pragma omp for 
        for (std::size_t i = 0; i < g_limit; ++i) { u[i] = rhs[i]; }
      }

      /* Retrieve timing at the beginning of Mx = g.                */
      #pragma omp single
      { 
        begin = timing::read();
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
        if constexpr (s == solver::dense_cholesky) {
          auto error = cholesky::solve(M_dense_cholesky[k], mgp, sbp);
          if (error) {
              std::cout << "Choleksy solve failed code " 
                  << *error << std::endl;
          }
        }
        else if constexpr (s == solver::sparse_qr) {
          status = mkl_sparse_d_qr_solve(
            SPARSE_OPERATION_NON_TRANSPOSE, M[td][k], nullptr,
            SPARSE_LAYOUT_COLUMN_MAJOR, 1, mgp , sbp.n * sbp.n, 
            gp, sbp.n * sbp.n);
          mkl_sparse_status(status);
        }

      }

      /* Calculate timing for Mx = g.                               */
      #pragma omp single
      {
        // __itt_task_end(domain);
        end = timing::read();
        trace.push_back(end - begin);
        logging::out << "Computed u = M^-1 * rhs." << std::endl;
      } 

      // #pragma omp single  
      // for (std::size_t i = 0; i < g_limit; ++i) {
      //   std::cout << Mg[i] << std::endl;
      // }

    } /* End OpenMP parallel region.                                */


    // logging::out << std::setw(14) << std::fixed << end - begin
    //  << " s # " << 
    logging::out << "Computed u (M \\ b)." << std::endl;
    
    // Free right-hand side memory. 
    mkl_free(rhs);
    
    // Gather solution on rank 0
    std::vector<int> recieve_counts (sbp.n_ranks);
    std::vector<int> displacement (sbp.n_ranks);
    for (std::size_t i = 0; i != sbp.n_ranks; ++i) {
      recieve_counts[i] = (i < alignment) 
        ? (sbp.n_blocks / sbp.n_ranks) + 1
        : (sbp.n_blocks / sbp.n_ranks);
      recieve_counts[i] *= sbp.n * sbp.n;
    }
    for (std::size_t i = 1; i < sbp.n_ranks; ++i) {
      displacement[i] = recieve_counts[i - 1] + displacement[i - 1];
    }

    real_t *uu = nullptr; 
    if(sbp.rank == 0) {
      sz = sizeof(real_t) * sbp.n * sbp.n * sbp.n_blocks; 
      uu = (double *) mkl_malloc(sz, 64);
    }
    begin = timing::read();

    sz = sbp.n * sbp.n * sbp.rank_limit_u;  
    MPI_Gatherv(
      &u[0],     // local send buffer 
      sz,         // # of elements to send
      MPI_DOUBLE, // type of element (vector, length varies per rank)
      uu,        // recieve buffer (only on rank 0)
      &recieve_counts[0], // list of # to recieve from each rank 
      &displacement[0],   // index of where to insert each vector
      MPI_DOUBLE,     // recieved type
      0,              // receiving rank 
      MPI_COMM_WORLD);

    // this is stored block by block but we need to print it global style 
    
    // TO PRINT U 
    
    /*
    for (std::size_t i = 0; i < sbp.n_blocks_dim; ++i) {
      std::size_t lim_a = (i == sbp.n_blocks_dim - 1) ? sbp.n : sbp.n - 1;
      for (std::size_t j = 0; j < lim_a; ++j) {
        for (std::size_t k = 0; k < sbp.n_blocks_dim; ++k) {
          std::size_t lim_b = (k == sbp.n_blocks_dim - 1) ? sbp.n : sbp.n - 1;
          for (std::size_t l = 0; l < lim_b; ++l) {
            auto index = 
                (i * sbp.n * sbp.n * sbp.n_blocks_dim) 
              + (k * sbp.n * sbp.n) 
              + (j * sbp.n) 
              + l;
            std::cout << uu[index] << " ";
          }
        }
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    */
    

    /*
    for (std::size_t i = 0; i < sbp.n; ++i) {
      for (std::size_t j = 0; j < sbp.n; ++j) {
        auto index = i * sbp.n + j;
        std::cout << uu[index] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
    for (std::size_t i = 0; i < sbp.n; ++i) {
      for (std::size_t j = 0; j < sbp.n; ++j) {
        auto index = (sbp.n * sbp.n) + (i * sbp.n + j);
        std::cout << uu[index] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    */

    //for (std::size_t i = 0; i< sbp.n * sbp.n * sbp.n_blocks; ++i) {
    //  std::cout << uu[i] << " ";
    //}
    //std::cout << std::endl;
    

    end = timing::read();
    trace.push_back(end - begin);
    // Cleanup everything we allocated.
    for (auto &e: M) {
      for (auto &ee: e)
        mkl_sparse_destroy(ee);
    }
    if (sbp.rank == 0) {
      mkl_free(Lambda_A);
    }
    mkl_sparse_destroy(D);
    // mkl_free(LAMBDA);
  
    
    mkl_free(u);

    mkl_free(boundary_solution);
    mkl_free(sources);
    for (auto &e: Fdense) mkl_free(e);
    for (auto &e: MF) mkl_free(e);
    mkl_free(Mg);
    
    if(sbp.rank == 0) {

        std::vector<std::string> key = {
          "LL^T=M            ",
          "F                 ",
          "X=M^-1F           ",
          "(Rank 0) LA=D-F^TX",
          "(Rank 0) LL^T=LA  ",
          "x=M^-1g           ",
          "Lb=-F^Tx          ",
          "Reduce sum Lb     ",
          "(Rank 0) L=L^-1Lb ",
          "Broadcast L       ",
          "b=g-FL            ",
          "u=M^-1b           ",
          "Gather u          "}; 


        std::size_t offset = trace.size();
        std::vector<real_t> traces (offset * sbp.n_ranks);
        MPI_Gather(&trace[0], offset, MPI_DOUBLE, &traces[0], offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // printf("Values collected on process %d: %d, %d, %d, %d.\n", my_rank, buffer[0], buffer[1], buffer[2], buffer[3]);
        
        double psolve = 0;
        double ssolve = 0;
        double pspmv  = 0;
        double pwork  = 0;
        std::cout << "                  ";
        for (std::size_t i = 0; i != sbp.n_ranks; ++i) { 
          std::cout << ",      Rank " << std::setw(2) << i;
        }
        
        for (std::size_t i = 0; i != offset; ++i) { 
          std::cout << std::endl << key[i];
          for (std::size_t j = 0; j != sbp.n_ranks; ++j) { 
            std::cout << ", " << std::setw(8) << std::scientific << traces[offset * j + i];
            if (i == 5 || i == 11) {
              psolve += traces[offset * j + i];
              pwork += traces[offset * j + i];
            }
            else if (i == 8) {
              ssolve += traces[offset * j + i];
            }
            else if (i == 10 || i == 6) {
              pspmv += traces[offset * j + i];
              pwork += traces[offset * j + i];
            }
          }
        }
        std::cout << std::endl;
        std::cout << std::endl;

        std::cout << "Parallel solve time : " 
          << std::setw(8) << std::scientific 
          << psolve << std::endl;
        std::cout << "Serial solve time   : " 
          << std::setw(8) << std::scientific 
          << ssolve << std::endl;
        std::cout << "Ratio               : " 
          << psolve / ssolve << std::endl;
        std::cout << "Ratio (Global)      : " 
          << pwork / ssolve << std::endl;
        std::cout << "Total running time  : " 
          << std::setw(8) << std::scientific 
          << pwork + ssolve << std::endl;
        std::cout << "                    : " 
          << std::fixed
          << pwork + ssolve << " sec" << std::endl;

        std::time_t ttt = std::time(0);  // t is an integer type

        std::ofstream file;
        std::stringstream filename;
        filename << "results/r_" << sbp.n_ranks << "_t_" << sbp.n_threads 
          << "_e_" << sbp.n_blocks << "_n2_" << sbp.n * sbp.n << "_"
          << ttt << ".csv";
        file.open(filename.str());
        file << "                ";
        for (std::size_t i = 0; i != sbp.n_ranks; ++i) { 
          file << ",      Rank " << std::setw(2) << i;
        }
        
        for (std::size_t i = 0; i != offset; ++i) { 
          file << std::endl << key[i];
          for (std::size_t j = 0; j != sbp.n_ranks; ++j) { 
            file << ", " << std::setw(8) << std::scientific << traces[offset * j + i];
          }
        }
        file << std::endl;
        file.close();
    }
    else {
      std::size_t offset = trace.size();
      MPI_Gather(&trace[0], offset, MPI_DOUBLE, NULL, offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }    
    
    return;
}

/* 

0.000676 s # Computed M.
      0.000374 s # Computed F.
      0.000319 s # Computed F.
      0.000291 s # Computed F.
      0.041461 s # Computed MX=F.
      0.044503 s # Computed MX=F.
      0.019742 s # Computed Lambda_A (D - FT * M \ F).
      0.000224 s # Computed Mx = g.
      0.000014 s # Computed LAMBDAb (gd - FT * M \ g).
      0.026061 s # Computed Lambda_A (D - FT * M \ F).
      0.000126 s # Computed Mx = g.
      0.000022 s # Computed LAMBDAb (gd - FT * M \ g).
      0.087782 s # Computed MX=F.
      0.026921 s # Factorized Lambda_A.
      0.027033 s # Factorized Lambda_A.
      0.000075 s # Computed LAMBDA (Lambda_A \ LAMBDAb).
      0.011034 s # Computed LAMBDA (Lambda_A \ LAMBDAb).
      0.011084 s # Computed b (g - F * LAMBDA).
      0.010832 s # Computed u (M \ b).
      0.010960 s # Computed b (g - F * LAMBDA).
      0.002984 s # Computed u (M \ b).
      0.047900 s # Computed Lambda_A (D - FT * M \ F).
      0.009134 s # Computed Mx = g.
      0.027810 s # Computed LAMBDAb (gd - FT * M \ g).
      0.000729 s # Factorized Lambda_A.
      0.000028 s # Computed LAMBDA (Lambda_A \ LAMBDAb).
      0.013223 s # Computed b (g - F * LAMBDA).
      0.008986 s # Computed u (M \ b).
*/