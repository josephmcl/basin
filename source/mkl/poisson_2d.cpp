#include "poisson_2d.h"

void poisson_2d::problem(std::size_t vln, std::size_t eln) {

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

    sbp.n_blocks = n_blocks; // additional non sbp-sat info. but  
    sbp.n_interfaces = n_interfaces; // useful to have along. 
    sbp.n_blocks_dim = l_blocks;
    sbp.n_threads = static_cast<std::size_t>(omp_get_max_threads());


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
    make_m(&M[0], sbp, {1, 1, 2, 1});
    make_m(&M[1], sbp, {1, 1, 1, 1});
    make_m(&M[2], sbp, {1, 1, 1, 2});
    
    // compute_lambda_matrix();

    std::cout << "computing f.. \n";
    // Compute F components.
    auto Fsparse = std::vector<sparse_matrix_t>(4);
    auto Fdense = std::vector<real_t *>(4);
    compute_f(Fsparse, Fdense, sbp);

    std::cout << "computed f\n";
    // Compute solve of MX = F.
    std::vector<real_t *> MF;
    compute_mf(MF, M, Fdense, sbp);
    
    // Cleanup everything we allocated.
    mkl_free(boundary_solution);
    mkl_free(sources);
    mkl_free(g);
    for (auto &e: M) mkl_free(e);
    for (auto &e: Fsparse) mkl_free(e);
    for (auto &e: Fdense) mkl_free(e);

    return;
}