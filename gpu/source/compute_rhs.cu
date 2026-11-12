#include "compute_rhs.cuh"
#include "cuda_check.cuh"
#include "spmv.cuh"

#include <stdexcept>

// Port of compute_rhs in hipbasin/basin/gpu/source/compute_rhs.cpp (the
// per-block version, driven by F_symbols) with explicit device memory.
// Differences from the HIP version, matching MKL's rhs = g - F^T lambda:
//
//  * Multiplies by F^T (n^2 x n). The HIP version multiplied by F (n x n^2):
//    it read n^2 entries from a length-n lambda segment and wrote only the
//    first n entries of the block's rhs.
//  * Subtracts (alpha = -1). The HIP version added (alpha = +1).
//  * Starts from rhs = g with a device-to-device copy, instead of adding g on
//    the host afterwards.
//
// MKL's compute_rhs is passed only this rank's slice of lambda but indexes it
// by global interface; this version takes the full lambda instead.
void compute_rhs(
    real_t                        *d_rhs,
    real_t const                  *d_g,
    std::vector<device_csr> const &d_FT,
    real_t const                  *d_lambda,
    vv<std::size_t> const         &F_symbols,
    components const              &sbp) {

  const std::size_t n2 = sbp.n * sbp.n;

  if (d_FT.size() != 4) {
    throw std::invalid_argument("compute_rhs: expected 4 transposed face operators");
  }
  for (auto const &ft : d_FT) {
    if (static_cast<std::size_t>(ft.n) != n2 ||
        static_cast<std::size_t>(ft.m) != sbp.n) {
      throw std::invalid_argument("compute_rhs: F^T is not n^2 x n");
    }
  }

  // rhs = g
  CUDA_CHECK(cudaMemcpy(d_rhs, d_g, sizeof(real_t) * n2 * sbp.rank_limit_u,
                        cudaMemcpyDeviceToDevice));

  for (std::size_t j = 0; j != sbp.rank_limit_u; ++j) {
    const std::size_t j_global = sbp.rank_index_u[j];
    for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
      const std::size_t face = F_symbols[j_global][i];
      if (face == 0) continue;

      // rhs_j -= F[face]^T * lambda_i   (rhs by local index, symbols by global)
      scalar_csr(d_FT[face - 1], d_lambda + sbp.n * i, d_rhs + n2 * j,
                 -1.0, 1.0);
    }
  }
}
