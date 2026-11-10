#pragma once

#include <vector>

#include "components.h"
#include "definitions.h"
#include "device_csr.cuh"

// Right-hand side of the trace system, on the device:
//
//   lambda_b_i = - sum over this rank's blocks j adjacent to interface i of
//                  F[face(i, j)] * Mg_j
//
//   d_lambda_b  device, n * n_interfaces, indexed by GLOBAL interface.
//               Overwritten (zeroed first). With several ranks each rank
//               holds a partial sum: reduce across ranks afterwards
//               (MPI_Allreduce with MPI_SUM).
//   d_F         F0..F3 on the device (W, E, S, N), each n x n^2.
//   d_Mg        device, M^-1 g, n^2 * sbp.rank_limit_u, LOCAL block index.
//   FT_symbols  [interface][global block]: 0 if not adjacent, otherwise the
//               1-based face of that block on the interface (1 W, 2 E, 3 S, 4 N).
//
// Launches on the default stream and returns without synchronizing.
void compute_lambda_b(
    real_t                        *d_lambda_b,
    std::vector<device_csr> const &d_F,
    real_t const                  *d_Mg,
    vv<std::size_t> const         &FT_symbols,
    components const              &sbp);
