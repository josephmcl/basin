#pragma once

#include "definitions.h"
#include "device_csr.cuh"

// Sparse matrix-vector products on the device:
//
//     y = alpha * A * x + beta * y
//
// x and y are device pointers and must not overlap. beta == 0 overwrites y
// without reading it. Both launch on the default stream and return without
// synchronizing.

// One thread per row. Best when rows hold only a few nonzeros, which is the
// case for every operator compute_g uses (B: 0-3 per row, H: 1 per row).
void scalar_csr(device_csr const &A, real_t const *x, real_t *y,
                real_t alpha, real_t beta);

// A group of 2-32 threads per row, sized from the average nonzeros per row,
// with a warp-shuffle reduction. Meant for denser rows.
void vector_csr(device_csr const &A, real_t const *x, real_t *y,
                real_t alpha, real_t beta);
