#include "poisson_2d.h"
#include "timing.h"

void poisson_2d::problem(std::size_t vln, std::size_t eln) {
    
    timing::init();

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

    mkl_set_num_threads(sbp.n_threads);

    MKL_INT *piv = (MKL_INT *) mkl_malloc(
      sizeof(MKL_INT) * sbp.n * sbp.n_interfaces, 64);
    // memset(piv, 0, sizeof(MKL_INT) * sbp.n * sbp.n_interfaces);
    for (std::size_t i = 0; i != sbp.n * sbp.n_interfaces; ++i) {
      piv[i] = 0;
    }
    double *λA = (double *) mkl_malloc(
      sizeof(double) * sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces, 64);
    for (std::size_t i = 0; i != sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces; ++i) {
      λA[i] = 0;
    }
    // memset(λA, 0, sizeof(double) * sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces);


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


    std::cout << "Square hybrid specs:" << std::endl 
    << " | local problem size: " << n_points_x << " x " << n_points_x 
    << std::endl << " | span: " << span 
    << std::endl << " | total blocks: " << l_blocks << " x " << l_blocks 
      << " (" << n_blocks << ") "
    << std::endl << " | total interfaces: " << n_interfaces 
    << std::endl << " | threads: " << sbp.n_threads 
    << std::endl;

    
    auto M = std::vector<sparse_matrix_t>(3);
    auto begin = timing::read();
    make_m(&M[0], sbp, {1, 1, 2, 1});
    make_m(&M[1], sbp, {1, 1, 1, 1});
    make_m(&M[2], sbp, {1, 1, 1, 2});
    auto end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin
      << " s # " << "Computed M." << std::endl;
    
    // compute_lambda_matrix();

    // Compute F components.
    auto Fsparse = std::vector<sparse_matrix_t>(4);
    auto Fdense = std::vector<real_t *>(4);
    begin = timing::read();
    compute_f(Fsparse, Fdense, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin 
    << " s # " << "Computed F." << std::endl;

    
    std::vector<real_t *> MF;
    MF.resize(M.size() * Fdense.size());
    for (std::size_t index = 0; index != MF.size(); ++index) {
        MF[index] = (real_t *) mkl_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.n, 64);
        memset(MF[index], 0, sizeof(real_t) * sbp.n * sbp.n * sbp.n);
    }
    
    // Compute solve of MX = F.
    begin = timing::read();
    compute_mf(MF, M, Fdense, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin 
    << " s # " << "Computed MX=F." << std::endl;

    // Setup D matrix.
    sparse_matrix_t D;
    compute_d(&D, sbp, interfaces);
    // std::cout << "computed d " << std::endl;

    // Setup interface list.
    vv<std::size_t> lambda_indices;
    make_interface_list(lambda_indices, F_symbols, FT_symbols, sbp);

    // Compute λA.
    begin = timing::read();
    // std::cout << sizeof(double) * sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces << std::endl;
    compute_lambda_a(λA, &D, Fsparse, MF, F_symbols, FT_symbols, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin 
      << " s # " << "Computed λA (D - FT * M \\ F)." << std::endl;

    /*
    sparse_index_base_t indexing;
    MKL_INT rows, cols, *rowst, *rowe, *coli, *ia;
    real_t *vals;
    auto status = mkl_sparse_d_export_csr(
      λA, &indexing, &rows, &cols, &rowst, &rowe, &coli, &vals);
    mkl_sparse_status(status);
    
    ia = (MKL_INT *) MKL_malloc(sizeof(MKL_INT) * ((sbp.n * sbp.n_interfaces) + 1), 64);
    std::memcpy(&ia[0], &rowst[0], sizeof(MKL_INT) * sbp.n * sbp.n_interfaces);
    ia[sbp.n * sbp.n_interfaces] = rowe[sbp.n * sbp.n_interfaces - 1];

    std::cout << rows << std::endl;
    std::cout << cols << std::endl;
    for (int i = 0; i < rowe[rows - 1]; ++i) {
      std::cout << vals[i] << " ";
    }
    std::cout << std::endl << std::endl;

    for (int i = 0; i < rowe[rows - 1]; ++i) {
      std::cout << coli[i] << " ";
    }
    std::cout << std::endl << std::endl;

    for (int i = 0; i < sbp.n * sbp.n_interfaces + 1; ++i) {
      std::cout << ia[i] << " ";
    }
    std::cout << std::endl;
    */

    // Compute Mx = g
    real_t *Mg;
    std::size_t sz = sizeof(real_t) * sbp.n * sbp.n * sbp.n_blocks;
    Mg = (double *) mkl_malloc(sz, 64);
    memset(Mg, 0, sz);

    begin = timing::read();
    compute_mg(Mg, M, g, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin 
    << " s # " << "Computed Mx = g." << std::endl;

    // Compute λb
    real_t *λb;
    sz = sizeof(real_t) * sbp.n * sbp.n_interfaces;
    λb = (double *) mkl_malloc(sz, 64);
    memset(λb, 0, sz);    

    begin = timing::read();
    compute_lambda_b(λb, Fsparse, Mg, FT_symbols, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin 
    << " s # " << "Computed λb (gd - FT * M \\ g)." << std::endl;

    // std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // std::cout << sizeof(MKL_INT) * sbp.n * sbp.n_interfaces << std::endl;
    // std::cout << sizeof(double) * sbp.n * sbp.n_interfaces * sbp.n * sbp.n_interfaces << std::endl;
    begin = timing::read();
    // Compute λ
    // std::cout << "potato" << std::endl;
    initialize_lambda(λA, piv, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin
      << " s # " << "Factorized λA." << std::endl;

    // real_t *λ;
    // sz = sizeof(real_t) * sbp.n * sbp.n_interfaces;
    // λ = (double *) mkl_malloc(sz, 64);
    // memset(λ, 0, sz);
    begin = timing::read();
    compute_lambda(λA, piv, λb, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin
      << " s # " << "Computed λ (λA \\ λb)." << std::endl;

    real_t *rhs;
    sz = sizeof(real_t) * sbp.n * sbp.n * sbp.n_blocks;
    rhs = (double *) mkl_malloc(sz, 64);
    memset(rhs, 0, sz);

    begin = timing::read();
    compute_rhs(rhs, Fsparse, λb, F_symbols, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin 
      << " s # " << "Computed b (g - F * λ)." << std::endl;

    real_t *u;
    sz = sizeof(real_t) * sbp.n * sbp.n * sbp.n_blocks;
    u = (double *) mkl_malloc(sz, 64);

    begin = timing::read();
    compute_u(u, M, rhs, sbp);
    end = timing::read();
    logging::out << std::setw(14) << std::fixed << end - begin
      << " s # " << "Computed u (M \\ b)." << std::endl;


    // Cleanup everything we allocated.
    for (auto &e: M) mkl_sparse_destroy(e);
    for (auto &e: Fsparse) mkl_sparse_destroy(e);
    mkl_free(λA);
    mkl_sparse_destroy(D);
    mkl_free(λb);
    // mkl_free(λ);
    mkl_free(rhs);
    mkl_free(u);

    mkl_free(boundary_solution);
    mkl_free(sources);
    mkl_free(g);
    for (auto &e: Fdense) mkl_free(e);
    for (auto &e: MF) mkl_free(e);
    mkl_free(Mg);
    

    return;
}