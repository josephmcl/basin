#include "compute_g.cuh"
#include "cuda_check.cuh"
#include "spmv.cuh"

#include <stdexcept>

// Port of hipbasin/basin/gpu/source/compute_g.cpp with explicit device
// memory. Differences from the HIP version:
//
//  * Boundary SpMVs accumulate (beta = 1) into a zeroed g. The HIP version
//    used beta = 0, so on a corner block (two boundary faces) the second
//    face's SpMV overwrote the first face's contribution.
//  * g is indexed by the LOCAL block and sources by the GLOBAL block. The HIP
//    version indexed g by the global block, which overruns the
//    rank_limit_u-sized buffer on any rank but 0.
//  * -H * source is one SpMV with alpha = -1 instead of negating into a
//    temporary first.
//  * No device sync after every launch: launches on one stream already run
//    in order.
void compute_g(
    real_t                        *d_g,
    std::vector<device_csr> const &d_boundaries,
    device_csr const              &d_h,
    real_t const                  *d_solutions,
    real_t const                  *d_sources,
    vv<std::size_t> const         &boundary_type_map,
    vv<std::size_t> const         &boundary_data_map,
    components const              &sbp) {

  const std::size_t n2 = sbp.n * sbp.n;
  const std::size_t face_size = sbp.n * sbp.n_blocks_dim;
  constexpr std::size_t faces = 4;  // blocks are rectangular

  if (d_boundaries.size() != 2 * faces) {
    throw std::invalid_argument("compute_g: expected 8 boundary operators");
  }
  for (auto const &b : d_boundaries) {
    if (static_cast<std::size_t>(b.n) != n2 ||
        static_cast<std::size_t>(b.m) != sbp.n) {
      throw std::invalid_argument("compute_g: boundary operator is not n^2 x n");
    }
  }
  if (static_cast<std::size_t>(d_h.n) != n2 ||
      static_cast<std::size_t>(d_h.m) != n2) {
    throw std::invalid_argument("compute_g: H is not n^2 x n^2");
  }

  CUDA_CHECK(cudaMemset(d_g, 0, sizeof(real_t) * n2 * sbp.rank_limit_u));

  for (std::size_t block = 0; block != sbp.rank_limit_u; ++block) {
    // block is the local index 0, 1, 2, ...; relblock is the global block
    // number k, k+1, k+2, ... used to look up maps and sources.
    const std::size_t relblock = sbp.rank_index_u[block];
    real_t *gi = d_g + n2 * block;

    for (std::size_t face = 0; face != faces; ++face) {
      // 0 = interior face (an interface), 1 = Dirichlet, 2 = Neumann.
      const std::size_t b_type = boundary_type_map[relblock][face];
      if (b_type == 0) continue;

      // Operator index: face + 4 * (type - 1).
      const std::size_t ti = face + (b_type - 1) * faces;

      // Which length-n segment of this face's edge belongs to this block:
      // [..face 0..][..face 1..][..face 2..][..face 3..], each n * L long.
      const std::size_t di = boundary_data_map[relblock][face] - 1;
      const real_t *solution = d_solutions + face * face_size + di * sbp.n;

      // gi += B[ti] * solution
      scalar_csr(d_boundaries[ti], solution, gi, 1.0, 1.0);
    }

    // gi -= H * source
    scalar_csr(d_h, d_sources + n2 * relblock, gi, -1.0, 1.0);
  }
}
