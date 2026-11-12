#include "compute_lambda_b.cuh"
#include "cuda_check.cuh"
#include "spmv.cuh"

#include <stdexcept>

// Port of compute_lambda_b in hipbasin/basin/gpu/source/compute_lambda_b.cpp
// (the per-block version, driven by FT_symbols) with explicit device memory.
// Same math; the only change is no device sync after every launch.
void compute_lambda_b(
    real_t                        *d_lambda_b,
    std::vector<device_csr> const &d_F,
    real_t const                  *d_Mg,
    vv<std::size_t> const         &FT_symbols,
    components const              &sbp) {

  const std::size_t n2 = sbp.n * sbp.n;

  if (d_F.size() != 4) {
    throw std::invalid_argument("compute_lambda_b: expected 4 face operators F0..F3");
  }
  for (auto const &f : d_F) {
    if (static_cast<std::size_t>(f.n) != sbp.n ||
        static_cast<std::size_t>(f.m) != n2) {
      throw std::invalid_argument("compute_lambda_b: F is not n x n^2");
    }
  }

  CUDA_CHECK(cudaMemset(d_lambda_b, 0, sizeof(real_t) * sbp.n * sbp.n_interfaces));

  // Two blocks touch each interface, so two launches accumulate into the same
  // lambda_b_i. That is safe only because launches on one stream run in order;
  // keep them on one stream (or use atomics) if this is ever parallelized.
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.rank_limit_u; ++j) {
      const std::size_t j_global = sbp.rank_index_u[j];
      const std::size_t face = FT_symbols[i][j_global];
      if (face == 0) continue;

      // lambda_b_i -= F[face] * Mg_j   (Mg by local index, symbols by global)
      scalar_csr(d_F[face - 1], d_Mg + n2 * j, d_lambda_b + sbp.n * i,
                 -1.0, 1.0);
    }
  }
}
