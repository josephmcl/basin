#pragma once

#include <vector>

#include "components.h"
#include "definitions.h"
#include "device_csr.cuh"

// Right-hand side g for every block this rank owns, computed on the device:
//
//   g_b = sum over b's physical-boundary faces f of  B[ti(f)] * g_f  -  H * source_b
//
//   d_g                device, n^2 * sbp.rank_limit_u, indexed by LOCAL block.
//                      Overwritten (zeroed first).
//   d_boundaries       B0..B7 on the device, each n^2 x n.
//                      [0..3] Dirichlet W,E,S,N; [4..7] Neumann W,E,S,N.
//   d_h                sbp.hl (2D quadrature H) on the device, n^2 x n^2.
//   d_solutions        device, boundary data laid out [W | E | S | N], each
//                      edge n * n_blocks_dim long.
//   d_sources          device, n^2 * n_blocks: ALL blocks, indexed by GLOBAL
//                      block number.
//   boundary_type_map  [block][face]: 0 interior, 1 Dirichlet, 2 Neumann.
//   boundary_data_map  [block][face]: 1-based segment index along that edge.
//
// Launches on the default stream and returns without synchronizing.
void compute_g(
    real_t                        *d_g,
    std::vector<device_csr> const &d_boundaries,
    device_csr const              &d_h,
    real_t const                  *d_solutions,
    real_t const                  *d_sources,
    vv<std::size_t> const         &boundary_type_map,
    vv<std::size_t> const         &boundary_data_map,
    components const              &sbp);
