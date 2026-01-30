#pragma once

// MPI rank distribution and setup utilities.

#include <mpi.h>
#include <vector>
#include <iostream>
#include <algorithm>

#include "components.h"

namespace mpi_setup {

// Initialize MPI parameters in the components struct.
// Distributes blocks and interfaces across MPI ranks.
inline void init(components &sbp) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    sbp.n_ranks = static_cast<std::size_t>(world_size);
    sbp.rank = static_cast<std::size_t>(world_rank);

    // Block distribution: alignment is the number of ranks with an extra block.
    std::size_t alignment = sbp.n_blocks % sbp.n_ranks;
    sbp.rank_limit_u = (sbp.rank < alignment)
        ? (sbp.n_blocks / sbp.n_ranks) + 1
        : (sbp.n_blocks / sbp.n_ranks);

    // Interface distribution
    std::size_t interface_alignment = sbp.n_interfaces % sbp.n_ranks;
    sbp.rank_limit_iu = (sbp.rank < interface_alignment)
        ? (sbp.n_interfaces / sbp.n_ranks) + 1
        : (sbp.n_interfaces / sbp.n_ranks);

    // Compute block indices for this rank
    std::size_t slack = (sbp.rank < alignment) ? 0 : std::min(sbp.rank, alignment);
    sbp.rank_index_u.clear();
    for (std::size_t i = 0; i != sbp.rank_limit_u; ++i) {
        sbp.rank_index_u.push_back(slack + sbp.rank_limit_u * sbp.rank + i);
    }

    // Compute interface indices for this rank
    slack = (sbp.rank < interface_alignment) ? 0 : std::min(sbp.rank, interface_alignment);
    sbp.rank_index_iu.clear();
    for (std::size_t i = 0; i != sbp.rank_limit_iu; ++i) {
        sbp.rank_index_iu.push_back(slack + sbp.rank_limit_iu * sbp.rank + i);
    }
}

// Print distribution info (rank 0 only).
inline void print_distribution(const components &sbp) {
    if (sbp.rank != 0) return;

    std::size_t alignment = sbp.n_blocks % sbp.n_ranks;

    std::cout << "Distributed solution range; blocks = "
              << sbp.n_blocks << "\n[";

    std::size_t cat = 0;
    for (std::size_t i = 0; i < sbp.n_ranks; ++i) {
        std::size_t r_limit, begin_range, end_range;
        std::size_t slack = (i < alignment) ? 0 : std::min(i, alignment);

        if (i < alignment) {
            r_limit = (sbp.n_blocks / sbp.n_ranks) + 1;
            begin_range = cat + slack + r_limit * i;
            end_range = begin_range + r_limit - 1;
            cat += 1;
        } else {
            r_limit = sbp.n_blocks / sbp.n_ranks;
            begin_range = cat + slack + r_limit * i;
            end_range = begin_range + r_limit - 1;
        }
        std::cout << begin_range << " ... " << end_range;
        if (i < sbp.n_ranks - 1) std::cout << " | ";
    }
    std::cout << "]\n";

    std::size_t interface_alignment = sbp.n_interfaces % sbp.n_ranks;
    std::cout << "Distributed interfaces range; interfaces = "
              << sbp.n_interfaces << "\n[";

    cat = 0;
    for (std::size_t i = 0; i < sbp.n_ranks; ++i) {
        std::size_t r_limit, begin_range, end_range;
        std::size_t slack = (i < interface_alignment) ? 0 : std::min(i, interface_alignment);

        if (i < interface_alignment) {
            r_limit = (sbp.n_interfaces / sbp.n_ranks) + 1;
            begin_range = cat + slack + r_limit * i;
            end_range = begin_range + r_limit - 1;
            cat += 1;
        } else {
            r_limit = sbp.n_interfaces / sbp.n_ranks;
            begin_range = cat + slack + r_limit * i;
            end_range = begin_range + r_limit - 1;
        }
        std::cout << begin_range << " ... " << end_range;
        if (i < sbp.n_ranks - 1) std::cout << " | ";
    }
    std::cout << "]\n";
}

// Compute receive counts and displacements for MPI_Gatherv.
inline void gather_params(const components &sbp,
                          std::vector<int> &recv_counts,
                          std::vector<int> &displacements) {
    std::size_t alignment = sbp.n_blocks % sbp.n_ranks;

    recv_counts.resize(sbp.n_ranks);
    displacements.resize(sbp.n_ranks);

    for (std::size_t i = 0; i != sbp.n_ranks; ++i) {
        recv_counts[i] = (i < alignment)
            ? (sbp.n_blocks / sbp.n_ranks) + 1
            : (sbp.n_blocks / sbp.n_ranks);
        recv_counts[i] *= sbp.n * sbp.n;
    }

    displacements[0] = 0;
    for (std::size_t i = 1; i < sbp.n_ranks; ++i) {
        displacements[i] = recv_counts[i - 1] + displacements[i - 1];
    }
}

} // namespace mpi_setup
