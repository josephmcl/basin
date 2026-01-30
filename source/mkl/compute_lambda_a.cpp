#include "compute_lambda_a.h"

#include <iostream>
#include <cmath>

namespace {

// Shared implementation for computing LambdaA.
// When distributed=true, uses rank-relative indexing for MPI distribution.
// When distributed=false, computes the full matrix on a single rank.
void compute_lambda_a_impl(
    double *lambdaA,
    sparse_matrix_t *D,
    std::vector<sparse_matrix_t> &F,
    std::vector<real_t*> &MF,
    vv<std::size_t> const &F_symbols,
    vv<std::size_t> const &FT_symbols,
    components &sbp,
    bool distributed) {

    matrix_descr da;
    da.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Allocate workspace for F^T * M^{-1} * F products
    std::size_t sz = sizeof(double) * sbp.n * sbp.n * F.size() * MF.size();
    double *FTMF = (double*)mkl_malloc(sz, 64);
    std::memset(FTMF, 0, sz);

    // Compute all F^T * M^{-1} * F block products
    for (std::size_t i = 0; i != F.size(); ++i) {
        for (std::size_t j = 0; j != MF.size(); ++j) {
            double *ftmf = &FTMF[(i * MF.size() + j) * sbp.n * sbp.n];

            sparse_status_t status = mkl_sparse_d_mm(
                SPARSE_OPERATION_NON_TRANSPOSE, 1., F[i], da,
                SPARSE_LAYOUT_COLUMN_MAJOR, MF[j],
                sbp.n, sbp.n * sbp.n,
                1., ftmf, sbp.n);
            mkl_sparse_status(status);
        }
    }

    // Build interface list: which (i,j,k) combinations contribute to LambdaA
    std::vector<std::size_t> interface_list;
    std::size_t slu = 0;

    // Determine iteration range based on distributed mode
    auto j_begin = distributed ? sbp.rank_index_iu.begin() : sbp.rank_index_iu.end();
    auto j_end = distributed ? sbp.rank_index_iu.end() : sbp.rank_index_iu.end();

    for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
        if (distributed) {
            // Distributed: iterate over this rank's interfaces
            for (auto j : sbp.rank_index_iu) {
                bool dirty = false;
                for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
                    if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {
                        if (!dirty) {
                            slu += 1;
                            dirty = true;
                        }
                        interface_list.push_back(i);
                        interface_list.push_back(j);
                        interface_list.push_back(k);
                        interface_list.push_back(slu - 1);
                    }
                }
            }
        } else {
            // Single rank: iterate over all interfaces
            for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
                bool dirty = false;
                for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
                    if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {
                        if (!dirty) {
                            slu += 1;
                            dirty = true;
                        }
                        interface_list.push_back(i);
                        interface_list.push_back(j);
                        interface_list.push_back(k);
                        interface_list.push_back(slu - 1);
                    }
                }
            }
        }
    }

    // Copy each block into the dense LambdaA matrix
    for (std::size_t i = 0; i < interface_list.size(); i += 4) {
        std::size_t a = interface_list[i];
        std::size_t b = interface_list[i + 1];
        std::size_t c = interface_list[i + 2];

        std::size_t ftindex = FT_symbols[a][c] - 1;
        std::size_t findex = F_symbols[c][b] - 1;

        // Determine M variant based on block position
        std::size_t r = (c % sbp.n_blocks_dim == 0) ? 0
                      : (c % sbp.n_blocks_dim == sbp.n_blocks_dim - 1) ? 2
                      : 1;

        std::size_t index = ((ftindex * MF.size()) + (r * 4) + findex) * sbp.n * sbp.n;
        double *ftmf = &FTMF[index];

        // Compute row offset for this rank
        std::size_t row_offset = distributed ? (b - sbp.rank_index_iu[0]) : b;

        for (std::size_t j = 0; j != sbp.n; ++j) {
            double *lat = &lambdaA[((j + (row_offset * sbp.n)) * sbp.n_interfaces * sbp.n) + (a * sbp.n)];
            for (std::size_t k = 0; k != sbp.n; ++k) {
                lat[k] -= ftmf[j * sbp.n + k];
            }
        }
    }

    // Add diagonal H matrix contribution
    std::size_t loop_limit = distributed ? sbp.rank_limit_iu * sbp.n : sbp.n_interfaces * sbp.n;
    std::size_t h_offset = distributed ? sbp.rank_index_iu[0] * sbp.n : 0;

    for (std::size_t j = 0; j != loop_limit; ++j) {
        std::size_t vector_index = j % sbp.n;
        std::size_t dense_index = (j * sbp.n_interfaces * sbp.n) + j + h_offset;
        lambdaA[dense_index] += sbp.h1v[vector_index] * 2.0 * sbp.TAU_VALUE;
    }

    mkl_free(FTMF);
}

} // anonymous namespace


// Distributed MPI version: each rank computes its portion of LambdaA.
void compute_lambda_a_mpi(
    double *lambdaA,
    sparse_matrix_t *D,
    std::vector<sparse_matrix_t> &F,
    std::vector<real_t*> &MF,
    vv<std::size_t> const &F_symbols,
    vv<std::size_t> const &FT_symbols,
    components &sbp) {

    compute_lambda_a_impl(lambdaA, D, F, MF, F_symbols, FT_symbols, sbp, true);
}


// Single-rank version: computes full LambdaA on one process.
void compute_lambda_a(
    double *lambdaA,
    sparse_matrix_t *D,
    std::vector<sparse_matrix_t> &F,
    std::vector<real_t*> &MF,
    vv<std::size_t> const &F_symbols,
    vv<std::size_t> const &FT_symbols,
    components &sbp) {

    compute_lambda_a_impl(lambdaA, D, F, MF, F_symbols, FT_symbols, sbp, false);
}
