

#include "poisson_2d_single.h"
//#include "timing.h"

void poisson_2d::single(components &sbp) {
    
    auto gw = [](real_t x, real_t y){return std::sin(π * x + π * y);};
    auto ge = [](real_t x, real_t y){return std::sin(π * x + π * y);};
    auto gs = [](real_t x, real_t y){return -π * std::cos(π * x + π * y);};
    auto gn = [](real_t x, real_t y){return π * std::cos(π * x + π * y);};

    // Multiply by 2 here because u_xx == u_yy
    auto source_function = [](real_t x, real_t y) {
        return -2. * π * π * sin(π * x + π * y);};

    auto exact_solution = [](real_t x, real_t y) {
        return sin(π * x + π * y);};

    std::vector<sparse_matrix_t> B;
    compute_b(B, sbp); 

    vv<std::size_t> boundary_data_map;
    vv<std::size_t> boundary_order_map;
    make_boundary_maps(boundary_data_map, boundary_order_map, sbp.n_blocks_dim);

    for (auto r : boundary_data_map) {
        for (auto e : r) {
        std::cout << std::setw(4);
        if (e == 0) std::cout << "   ·";
        else std::cout << e;
        }
        std::cout << std::endl;
    }

    for (auto r : boundary_order_map) {
        for (auto e : r) {
        std::cout << std::setw(4);
        if (e == 0) std::cout << "   ·";
        else std::cout << e;
        }
        std::cout << std::endl;
    }

    // Generate ranges for the x and y of each block, given the number of 
    // blocks in each dimension.
    auto block_grid = range_t(0, 1., sbp.n_blocks_dim + 1);
    std::vector<range_t> grids; 
    for (std::size_t i = 0; i != block_grid.size() - 1; ++i)  {
        grids.push_back(range_t(*block_grid[i], *block_grid[i + 1], sbp.
        n));
    }

    // Generate the solution at the boundary. This implementation generates 
    // vectors from range functions. 
    real_t *boundary_solution;
    compute_boundary_solution(  
        &boundary_solution, grids, {gw, ge, gs, gn}, {0., 1., 0., 1.});

    real_t *sources;
    compute_sources(&sources, grids, source_function);

    for (std::size_t i = 0; i != sbp.n; ++i) {
        for (std::size_t j = 0; j != sbp.n; ++j) {
            std::cout << source_function(static_cast<real_t>(i) / static_cast<real_t>(sbp.n - 1), static_cast<real_t>(j) / static_cast<real_t>(sbp.n - 1)) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "g" << std::endl;
    for (std::size_t i = 0; i != sbp.n; ++i) {
        std::cout << gw(0, static_cast<real_t>(i) / static_cast<real_t>(sbp.n - 1)) << " ";
    }
    std::cout << std::endl;
    for (std::size_t i = 0; i != sbp.n; ++i) {
        std::cout << ge(1, static_cast<real_t>(i) / static_cast<real_t>(sbp.n - 1)) << " ";
    }
    std::cout << std::endl;
    for (std::size_t i = 0; i != sbp.n; ++i) {
        std::cout << gs(static_cast<real_t>(i) / static_cast<real_t>(sbp.n - 1), 0) << " ";
    }
    std::cout << std::endl;
    for (std::size_t i = 0; i != sbp.n; ++i) {
        std::cout << gn(static_cast<real_t>(i) / static_cast<real_t>(sbp.n - 1), 1) << " ";
    }
    std::cout << std::endl;

    for (std::size_t i = 0; i != sbp.n; ++i) {
        for (std::size_t j = 0; j != sbp.n; ++j) {
            std::cout << sources[i * sbp.n + j] << " ";
        }
        std::cout << std::endl;
    }
    

    // Compute the hybrid system g terms. 
    real_t *g;
    compute_g(&g, B, boundary_solution, sources, boundary_order_map, 
    boundary_data_map, sbp); 

    for (std::size_t i = 0; i != sbp.n; ++i) {
        for (std::size_t j = 0; j != sbp.n; ++j) {
            std::cout << g[i * sbp.n + j] << " ";
        }
        // std::cout << std::endl;
    }
    std::cout << std::endl;

    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < sbp.n; ++j) {
            std::cout << boundary_solution[i * 5 + j] << " ";
        }
        std::cout << std::endl;
    }
     std::cout << std::endl;

    sparse_matrix_t A;
    make_m(&A, sbp, {1, 1, 2, 2});

    sparse_status_t status;
    matrix_descr dc;
    dc.type = SPARSE_MATRIX_TYPE_GENERAL;

    status = mkl_sparse_qr_reorder(A, dc);
    mkl_sparse_status(status);

    status = mkl_sparse_d_qr_factorize(A, nullptr);
    mkl_sparse_status(status);

    real_t *u = (real_t *) mkl_malloc(
          sizeof(double) * sbp.n * sbp.n, 64);
    for (std::size_t i = 0; i != sbp.n * sbp.n; ++i) {
        u[i] = 0;
    }
    
    status = mkl_sparse_d_qr_solve(
        SPARSE_OPERATION_NON_TRANSPOSE, A, nullptr,
        SPARSE_LAYOUT_COLUMN_MAJOR, 1, u , sbp.n * sbp.n, 
        g, sbp.n * sbp.n);
    mkl_sparse_status(status);

    for (std::size_t i = 0; i != sbp.n; ++i) {
        for (std::size_t j = 0; j != sbp.n; ++j) {
            std::cout << u[i * sbp.n + j] << " ";
        }
        // std::cout << std::endl;
    }
    std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;

    /*
    double *errv = (double *) malloc(sizeof(double) * sbp.n * sbp.n);
    double *err = (double *) malloc(sizeof(double) * sbp.n * sbp.n);
    for (std::size_t i = 0; i != sbp.n; ++i) {
        for (std::size_t j = 0; j != sbp.n; ++j) {
            errv[i * sbp.n + j] = u[i * sbp.n + j] - exact_solution(static_cast<real_t>(i) / static_cast<real_t>(sbp.n - 1), static_cast<real_t>(j) / static_cast<real_t>(sbp.n - 1));
        }
    }
    
    sparse_matrix_t h;
    sbp.hl.mkl(&h);

    matrix_descr da;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;
    status = mkl_sparse_d_mv(
        SPARSE_OPERATION_NON_TRANSPOSE, 1., 
        h, da,
        errv, 
        1., err);
    
    mkl_sparse_status(status);

    for (std::size_t i = 0; i != sbp.n; ++i) {
        for (std::size_t j = 0; j != sbp.n; ++j) {
            std::cout << err[i * sbp.n + j] << " ";
        }
    }
    std::cout << std::endl << std::endl;

    double errr;
    errr = cblas_ddot(sbp.n * sbp.n, err, 1, err, 1);

    std::cout << "error " << errr << std::endl;
    */
    std::cout << std::endl;
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