#include "poisson_2d.h"
#include "timing.h"

void poisson_2d::problem(std::size_t vln, std::size_t eln) {
    
    timing::init();

    // omp_set_nested(1);

    const std::size_t l_blocks = eln;
    const std::size_t n_blocks = l_blocks * l_blocks;
    auto span = 1. / static_cast<double>(l_blocks);
    auto n_points_x = vln;
    auto n = vln;

    auto sbp = components{n, span};

    sbp.τ = 42.; // hard code these coeffs for now. 
    sbp.β = 1.;

    auto gw = [](real_t x, real_t y){return std::sin(π * x + π * y);};
    auto ge = [](real_t x, real_t y){return std::sin(π * x + π * y);};
    auto gs = [](real_t x, real_t y){return -π * std::cos(π * x + π * y);};
    auto gn = [](real_t x, real_t y){return π * std::cos(π * x + π * y);};

    // Multiply by 2 here because u_xx == u_yy
    auto source_function = [](real_t x, real_t y) {
        return 2. * (-π * π * sin(π * x + π * y));};

    vv<std::size_t> interfaces;
    std::size_t n_interfaces = make_connectivity(interfaces, l_blocks);

    /*
    std::cout << "   Connectivity (nnz = " << n_interfaces 
              << ")" << std::endl;
    for (auto r : interfaces) {
        for (auto e : r) {
        std::cout << std::setw(4);
        if (e == 0) std::cout << "   ·";
        else std::cout << e;
        }
        std::cout << std::endl;
    }
    */

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
    std::size_t rank_limit_max = (sbp.n_blocks / sbp.n_ranks) + 1;

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

    // Compute flux components used by b. 
    std::vector<sparse_matrix_t> B;
    compute_b(B, sbp);  

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
        &boundary_solution, grids, {gw, ge, gs, gn}, {0., 1., 0., 1.});

    real_t *sources;
    compute_sources(&sources, grids, source_function);

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
    M.resize(sbp.n_threads);
    for (auto &e: M) {
      e.resize(3);
    }
    auto begin = timing::read();
    make_m(&M[0][0], sbp, {1, 1, 2, 1});
    make_m(&M[0][1], sbp, {1, 1, 1, 1});
    make_m(&M[0][2], sbp, {1, 1, 1, 2});
    auto end = timing::read();
    trace.push_back(end - begin);
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
    for (std::size_t i = 0; i != M.size(); ++i) {
      for (std::size_t j = 0; j != M[i].size(); ++j) {
      
        status = mkl_sparse_qr_reorder(M[i][j], dc);
        mkl_sparse_status(status);

        status = mkl_sparse_d_qr_factorize(M[i][j], nullptr);
        mkl_sparse_status(status);
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
    
    // Compute solve of MX = F.
    begin = timing::read();
    compute_mf(MF, M, Fdense, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin 
    // << " s # " << "Computed MX=F." << std::endl;

    // Setup D matrix.
    sparse_matrix_t D;
    compute_d(&D, sbp, interfaces);
    // std::cout << "computed d " << std::endl;

    // Setup interface list.
    vv<std::size_t> lambda_indices;
    make_interface_list(lambda_indices, F_symbols, FT_symbols, sbp);

    // Allocate memory for λA.
    double *λA = (double *) mkl_malloc(
      sizeof(double) * sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces, 64);
    for (std::size_t i = 0; i != sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces; ++i) {
      λA[i] = 0;
    }

    // Compute λA.
    begin = timing::read();
    // std::cout << sizeof(double) * sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces << std::endl;
    compute_lambda_a(λA, &D, Fsparse, MF, F_symbols, FT_symbols, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin 
    //  << " s # " << "Computed λA (D - FT * M \\ F)." << std::endl;

    // Compute Mx = g
    real_t *Mg;
    std::size_t sz = sizeof(real_t) * sbp.n * sbp.n * sbp.n_blocks;
    Mg = (double *) mkl_malloc(sz, 64);
    memset(Mg, 0, sz);

    begin = timing::read();
    compute_mg(Mg, M, g, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin 
    // << " s # " << "Computed Mx = g." << std::endl;

    
    real_t *λb;
    sz = sizeof(real_t) * sbp.n * sbp.n_interfaces;
    λb = (double *) mkl_malloc(sz, 64);
    memset(λb, 0, sz);    

    // Compute λb
    begin = timing::read();
    compute_lambda_b(λb, Fsparse, Mg, FT_symbols, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin 
    // << " s # " << "Computed λb (gd - FT * M \\ g)." << std::endl;

    // Allocate pivot table memory.
    MKL_INT *piv = (MKL_INT *) mkl_malloc(
      sizeof(MKL_INT) * sbp.n * sbp.n_interfaces, 64);
    // memset(piv, 0, sizeof(MKL_INT) * sbp.n * sbp.n_interfaces);
    for (std::size_t i = 0; i != sbp.n * sbp.n_interfaces; ++i) {
      piv[i] = 0;
    }

    // Compute λ
    begin = timing::read();
    initialize_lambda(λA, piv, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin
    //  << " s # " << "Factorized λA." << std::endl;

    begin = timing::read();
    compute_lambda(λA, piv, λb, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin
    //  << " s # " << "Computed λ (λA \\ λb)." << std::endl;

    // NOTE NOTE NOTE NOTE

    // Allocated right-hand side memory. 
    real_t *rhs;
    sz = sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u;
    rhs = (double *) mkl_malloc(sz, 64);
    memset(rhs, 0, sz);

    // Compute right-hand side. b = g - F * λ.
    begin = timing::read();
    compute_rhs(rhs, g, Fsparse, λb, F_symbols, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin 
    //   << " s # " << "Computed b (g - F * λ)." << std::endl;

    // Free g, F (sparse), and λ.
    mkl_free(g);
    for (auto &e: Fsparse) 
      mkl_sparse_destroy(e);
    mkl_free(λb);

    // Allocate solution memory. 
    real_t *u; 
    // NOTE: rank_limit_max here so we can gather the same size on 
    //       everything.
    sz = sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u; 
    u = (double *) mkl_malloc(sz, 64);
    memset(u, 0, sz);

    // Compute solution. u = M \ b.
    begin = timing::read();
    compute_u(u, M, rhs, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    // logging::out << std::setw(14) << std::fixed << end - begin
    //  << " s # " << "Computed u (M \\ b)." << std::endl;

    // Free right-hand side memory. 
    mkl_free(rhs);
    begin = timing::read();
    compute_u(u, M, rhs, sbp);
    end = timing::read();
    trace.push_back(end - begin);
    
    // Gather solution on rank 0
    /*
    std::vector<int> recieve_counts (sbp.n_ranks);
    std::vector<int> displacement (sbp.n_ranks);
    for (std::size_t i = 0; i != sbp.n_ranks; ++i) {
      recieve_counts[i] = (i < alignment) 
        ? (sbp.n_blocks / sbp.n_ranks) + 1
        : (sbp.n_blocks / sbp.n_ranks);
    }
    for (std::size_t i = 1; i < sbp.n_ranks; ++i) {
      displacement[i] = recieve_counts[i - 1] + displacement[i - 1];
    }

    real_t *uu = nullptr; 
    
    if(sbp.rank == 0) {
      begin = timing::read();
      sz = sizeof(real_t) * sbp.n * sbp.n * sbp.n_blocks; 
      uu = (double *) mkl_malloc(sz, 64);
    }

    sz = sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u; 
    MPI_Datatype send_type; 
    MPI_Type_vector(1, sz, 0, MPI_DOUBLE, &send_type);
    
    if(sbp.rank == 0) {
      MPI_Gatherv(
        &u[0],     // local send buffer 
        1,         // # of elements to send
        send_type, // type of element (vector, length varies per rank)
        uu,        // recieve buffer (only on rank 0)
        &recieve_counts[0], // list of # to recieve from each rank 
        &displacement[0],   // index of where to insert each vector
        MPI_DOUBLE,     // recieved type
        0,              // receiving rank 
        MPI_COMM_WORLD);
    }
    else {
      MPI_Gatherv(
        &u[0],     // local send buffer 
        1,         // # of elements to send
        send_type, // type of element (vector, length varies per rank)
        NULL,        // recieve buffer (only on rank 0)
        &recieve_counts[0], // list of # to recieve from each rank 
        &displacement[0],   // index of where to insert each vector
        MPI_DOUBLE,     // recieved type
        0,              // receiving rank 
        MPI_COMM_WORLD);
    }

    if(sbp.rank == 0) {
      end = timing::read();
      std::cout << "Gather took " << end - begin << std::endl;
    }
    */
    
    std::cout << "banana" << std::endl;



    // Cleanup everything we allocated.
    for (auto &e: M) {
      for (auto &ee: e)
        mkl_sparse_destroy(ee);
    }
    mkl_free(λA);
    mkl_sparse_destroy(D);
    // mkl_free(λ);
    mkl_free(u);

    mkl_free(boundary_solution);
    mkl_free(sources);
    for (auto &e: Fdense) mkl_free(e);
    for (auto &e: MF) mkl_free(e);
    mkl_free(Mg);
    
    if(sbp.rank == 0) {

        std::vector<std::string> key = {
          "LL^T=M   ",
          "F        ",
          "X=M^-1F  ",
          "λA=D-F^TX",
          "x=M^-1g  ",
          "λb=-F^Tx ",
          "LL^T=λA  ",
          "λ=λ^-1λb ",
          "b=g-Fλ   ",
          "u=M^-1b  "};

        std::size_t offset = trace.size();
        std::vector<real_t> traces (10 * sbp.n_ranks);
        MPI_Gather(&trace[0], 10, MPI_DOUBLE, &traces[0], 10, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // printf("Values collected on process %d: %d, %d, %d, %d.\n", my_rank, buffer[0], buffer[1], buffer[2], buffer[3]);

        std::cout << "         ";
        for (std::size_t i = 0; i != sbp.n_ranks; ++i) { 
          std::cout << ",      Rank " << std::setw(2) << i;
        }
        
        for (std::size_t i = 0; i != offset; ++i) { 
          std::cout << std::endl << key[i];
          for (std::size_t j = 0; j != sbp.n_ranks; ++j) { 
            std::cout << ", " << std::setw(8) << std::scientific << traces[offset * j + i];
          }
        }
        std::cout << std::endl;
    }
    else {
        MPI_Gather(&trace[0], 10, MPI_DOUBLE, NULL, 10, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
      0.019742 s # Computed λA (D - FT * M \ F).
      0.000224 s # Computed Mx = g.
      0.000014 s # Computed λb (gd - FT * M \ g).
      0.026061 s # Computed λA (D - FT * M \ F).
      0.000126 s # Computed Mx = g.
      0.000022 s # Computed λb (gd - FT * M \ g).
      0.087782 s # Computed MX=F.
      0.026921 s # Factorized λA.
      0.027033 s # Factorized λA.
      0.000075 s # Computed λ (λA \ λb).
      0.011034 s # Computed λ (λA \ λb).
      0.011084 s # Computed b (g - F * λ).
      0.010832 s # Computed u (M \ b).
      0.010960 s # Computed b (g - F * λ).
      0.002984 s # Computed u (M \ b).
      0.047900 s # Computed λA (D - FT * M \ F).
      0.009134 s # Computed Mx = g.
      0.027810 s # Computed λb (gd - FT * M \ g).
      0.000729 s # Factorized λA.
      0.000028 s # Computed λ (λA \ λb).
      0.013223 s # Computed b (g - F * λ).
      0.008986 s # Computed u (M \ b).
*/