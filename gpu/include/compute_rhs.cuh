#pragma once

#include <vector>

#include "components.h"
#include "definitions.h"
#include "device_csr.cuh"

// Right-hand side of the local block solves, on the device:
//
//   rhs_j = g_j - sum over interfaces i adjacent to block j of
//                 F[face(j, i)]^T * lambda_i
//
//   d_rhs      device, n^2 * sbp.rank_limit_u, LOCAL block index. Overwritten.
//   d_g        device, n^2 * sbp.rank_limit_u, LOCAL block index.
//   d_FT       F0^T..F3^T on the device (W, E, S, N), each n^2 x n.
//              Build on the host with csr_t(F[k], true), then to_device().
//   d_lambda   device, the FULL trace solution, n * n_interfaces, indexed by
//              GLOBAL interface (every rank needs all of it).
//   F_symbols  [global block][interface]: 0 if not adjacent, otherwise the
//              1-based face of the block on that interface (1 W, 2 E, 3 S, 4 N).
//
// Launches on the default stream and returns without synchronizing.
void compute_rhs(
    real_t                        *d_rhs,
    real_t const                  *d_g,
    std::vector<device_csr> const &d_FT,
    real_t const                  *d_lambda,
    vv<std::size_t> const         &F_symbols,
    components const              &sbp);
