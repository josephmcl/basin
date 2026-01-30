// Entry point for generating and storing precomputed solver components.
// This computes the M matrices and LambdaA (Schur complement) and saves them
// to disk for later use, separating factorization from solve.
//
// Usage: mpirun -np <ranks> basin-generate <n_points> <n_blocks_dim> <output_dir>

#include "poisson_2d.h"
#include "csr.h"
#include "serialize.h"
#include "mpi_setup.h"
#include "compute_m.h"
#include "compute_b.h"
#include "compute_f.h"
#include "compute_f_symbols.h"
#include "compute_mf.h"
#include "compute_lambda_a.h"
#include "compute_d.h"
#include "connect.h"
#include "timing.h"

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

void generate_and_save(std::size_t n_points, std::size_t n_blocks_dim,
                       const std::string &output_dir);

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <n_points> <n_blocks_dim> <output_dir>\n";
        return 1;
    }

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    std::size_t n_points = std::stoul(argv[1]);
    std::size_t n_blocks_dim = std::stoul(argv[2]);
    std::string output_dir = argv[3];

    generate_and_save(n_points, n_blocks_dim, output_dir);

    MPI_Finalize();
    return 0;
}

void generate_and_save(std::size_t n_points, std::size_t n_blocks_dim,
                       const std::string &output_dir) {
    timing::init();
    omp_set_nested(true);
    omp_set_max_active_levels(4);

    const std::size_t n_blocks = n_blocks_dim * n_blocks_dim;
    const double span = 1.0 / static_cast<double>(n_blocks_dim);

    // Initialize SBP components
    auto sbp = components{n_points, span};

    // Compute penalty parameters
    double space = 0.5 * (span / (n_points - 1));
    sbp.TAU_VALUE = (2.0 / span) + (2.0 / (span * (space / span)));
    sbp.BETA_VALUE = 1.0;

    sbp.n_blocks = n_blocks;
    sbp.n_blocks_dim = n_blocks_dim;
    sbp.n_threads = static_cast<std::size_t>(omp_get_max_threads());

    // Compute connectivity
    vv<std::size_t> interfaces;
    std::size_t n_interfaces = make_connectivity(interfaces, n_blocks_dim);
    sbp.n_interfaces = n_interfaces;

    // Initialize MPI distribution
    mpi_setup::init(sbp);
    mpi_setup::print_distribution(sbp);

    if (sbp.rank == 0) {
        std::cout << "\nGenerating components:\n"
                  << "  n_points: " << n_points << "\n"
                  << "  n_blocks: " << n_blocks << " (" << n_blocks_dim << "x"
                  << n_blocks_dim << ")\n"
                  << "  n_interfaces: " << n_interfaces << "\n"
                  << "  tau: " << sbp.TAU_VALUE << "\n"
                  << "  beta: " << sbp.BETA_VALUE << "\n"
                  << "  output: " << output_dir << "\n\n";

        // Create output directory
        mkdir(output_dir.c_str(), 0755);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // =========================================================================
    // Compute M matrices (3 variants for boundary conditions)
    // =========================================================================
    if (sbp.rank == 0) std::cout << "Computing M matrices..." << std::flush;

    auto M = vv<sparse_matrix_t>();
    M.resize(sbp.n_threads);
    for (auto &e : M) e.resize(3);

    double begin = timing::read();
    make_m(&M[0][0], sbp, {1, 1, 2, 1});  // West Neumann
    make_m(&M[0][1], sbp, {1, 1, 1, 1});  // All Dirichlet
    make_m(&M[0][2], sbp, {1, 1, 1, 2});  // East Neumann
    double end = timing::read();

    if (sbp.rank == 0)
        std::cout << " done (" << end - begin << " s)\n";

    // Copy M matrices for each thread
    matrix_descr dc;
    dc.type = SPARSE_MATRIX_TYPE_GENERAL;
    for (std::size_t i = 1; i != M.size(); ++i) {
        mkl_sparse_copy(M[0][0], dc, &M[i][0]);
        mkl_sparse_copy(M[0][1], dc, &M[i][1]);
        mkl_sparse_copy(M[0][2], dc, &M[i][2]);
    }

    // Factor M matrices
    if (sbp.rank == 0) std::cout << "Factoring M matrices..." << std::flush;
    begin = timing::read();
    for (std::size_t i = 0; i != M.size(); ++i) {
        for (std::size_t j = 0; j != M[i].size(); ++j) {
            mkl_sparse_qr_reorder(M[i][j], dc);
            mkl_sparse_d_qr_factorize(M[i][j], nullptr);
        }
    }
    end = timing::read();
    if (sbp.rank == 0)
        std::cout << " done (" << end - begin << " s)\n";

    // =========================================================================
    // Compute F matrices
    // =========================================================================
    if (sbp.rank == 0) std::cout << "Computing F matrices..." << std::flush;
    begin = timing::read();

    std::vector<sparse_matrix_t> Fsparse(4);
    std::vector<real_t*> Fdense(4);
    compute_f(Fsparse, Fdense, sbp);

    end = timing::read();
    if (sbp.rank == 0)
        std::cout << " done (" << end - begin << " s)\n";

    // =========================================================================
    // Compute MF = M^{-1} F
    // =========================================================================
    if (sbp.rank == 0) std::cout << "Computing M^{-1}F..." << std::flush;
    begin = timing::read();

    std::vector<real_t*> MF;
    compute_mf(MF, M, Fdense, sbp);

    end = timing::read();
    if (sbp.rank == 0)
        std::cout << " done (" << end - begin << " s)\n";

    // =========================================================================
    // Compute F symbols (connectivity)
    // =========================================================================
    vv<std::size_t> F_symbols(n_blocks, std::vector<std::size_t>(n_interfaces, 0));
    vv<std::size_t> FT_symbols(n_interfaces, std::vector<std::size_t>(n_blocks, 0));
    compute_f_symbols(F_symbols, FT_symbols, interfaces, sbp);

    // =========================================================================
    // Compute D matrix
    // =========================================================================
    sparse_matrix_t D;
    compute_d(&D, sbp, interfaces);

    // =========================================================================
    // Compute LambdaA (Schur complement)
    // =========================================================================
    if (sbp.rank == 0) std::cout << "Computing LambdaA..." << std::flush;
    begin = timing::read();

    std::size_t lambda_size = sbp.n * sbp.n_interfaces;
    double *LAMBDAA = (double*)mkl_malloc(
        sizeof(double) * lambda_size * sbp.n * sbp.rank_limit_iu, 64);
    std::memset(LAMBDAA, 0,
        sizeof(double) * lambda_size * sbp.n * sbp.rank_limit_iu);

    compute_lambda_a_mpi(LAMBDAA, &D, Fsparse, MF, F_symbols, FT_symbols, sbp);

    end = timing::read();
    if (sbp.rank == 0)
        std::cout << " done (" << end - begin << " s)\n";

    // =========================================================================
    // Convert to sparse CSR and save
    // =========================================================================
    if (sbp.rank == 0) std::cout << "Converting to sparse format..." << std::flush;
    begin = timing::read();

    csr<double> lambda_a_csr(lambda_size, lambda_size);
    for (std::size_t i = 0; i != sbp.rank_limit_iu * sbp.n; ++i) {
        for (std::size_t j = 0; j != lambda_size; ++j) {
            double val = LAMBDAA[i * lambda_size + j];
            if (val != 0) {
                lambda_a_csr(val, i + (sbp.rank_index_iu[0] * sbp.n), j);
            }
        }
    }

    end = timing::read();
    if (sbp.rank == 0)
        std::cout << " done (" << end - begin << " s)\n";

    // =========================================================================
    // Save to disk
    // =========================================================================
    if (sbp.rank == 0) {
        std::cout << "Saving to disk..." << std::flush;
        begin = timing::read();

        // Save parameters
        serialize::write_params(output_dir + "/params.txt",
                                n_points, n_blocks_dim,
                                sbp.TAU_VALUE, sbp.BETA_VALUE);

        // Save LambdaA (rank 0 gathers all pieces)
        // For now, each rank saves its portion
        end = timing::read();
        std::cout << " done (" << end - begin << " s)\n";
    }

    // Each rank saves its portion of LambdaA
    std::stringstream filename;
    filename << output_dir << "/lambda_a_rank" << sbp.rank << ".bin";
    serialize::write_csr(filename.str(), lambda_a_csr);

    if (sbp.rank == 0) {
        std::cout << "\nComponents saved to " << output_dir << "/\n"
                  << "  params.txt - Problem parameters\n"
                  << "  lambda_a_rank*.bin - LambdaA sparse matrix (per rank)\n";
    }

    // Cleanup
    mkl_free(LAMBDAA);
    for (auto &e : Fsparse) mkl_sparse_destroy(e);
    for (auto &e : Fdense) mkl_free(e);
    for (auto &e : MF) mkl_free(e);
    for (auto &e : M) {
        for (auto &ee : e) mkl_sparse_destroy(ee);
    }
    mkl_sparse_destroy(D);
}
