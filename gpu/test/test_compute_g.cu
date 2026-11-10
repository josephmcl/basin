// Checks the GPU compute_g against a straightforward host implementation of
// the same formula:
//
//   g_b = sum over b's boundary faces f of  B[ti(f)] * g_f  -  H * source_b
//
// B0..B7 come from ../operator/V_*_N_<n>_L_<L>/B<i>.* when the MKL generator
// has exported them (with .ndim/.mdim); otherwise random sparse n^2 x n
// matrices are used, which still exercises the kernels, the transfers and the
// face/segment/block indexing.
//
// Every case is repeated for 1, 2 and 3 simulated MPI ranks so that the
// local-vs-global block indexing is checked too.
//
// Usage: ./bin/test_compute_g [n=8] [L=3]

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include "components.h"
#include "compute_boundary_solution.h"
#include "compute_g.cuh"
#include "compute_sources.h"
#include "connect.h"
#include "cuda_check.cuh"
#include "data.h"
#include "definitions.h"
#include "device_csr.cuh"

namespace {

// Same block distribution as main()/poisson_2d::problem.
void assign_blocks(components &sbp, std::size_t n_ranks, std::size_t rank) {
  sbp.n_ranks = n_ranks;
  sbp.rank = rank;
  const std::size_t alignment = sbp.n_blocks % n_ranks;
  sbp.rank_limit_u = (rank < alignment) ? sbp.n_blocks / n_ranks + 1
                                        : sbp.n_blocks / n_ranks;
  const std::size_t slack = (rank < alignment) ? 0 : std::min(rank, alignment);
  sbp.rank_index_u.clear();
  for (std::size_t i = 0; i != sbp.rank_limit_u; ++i)
    sbp.rank_index_u.push_back(slack + sbp.rank_limit_u * rank + i);
}

// y += alpha * A * x on the host.
void spmv_host(csr_t const &A, real_t const *x, real_t *y, real_t alpha) {
  for (std::size_t i = 0; i < A.n; ++i) {
    real_t sum = 0;
    for (int j = A.r[i]; j < A.r[i + 1]; ++j) sum += A.v[j] * x[A.c[j]];
    y[i] += alpha * sum;
  }
}

// Random n^2 x n matrix: ~30% of rows get 1-3 nonzeros in distinct columns.
csr_t random_boundary_operator(std::size_t n, std::mt19937_64 &rng) {
  std::uniform_real_distribution<real_t> value(0.5, 1.5);
  std::uniform_int_distribution<int> coin(0, 9), count(1, 3), sign(0, 1);
  std::uniform_int_distribution<std::size_t> column(0, n - 1);
  csr_t B(n * n, n);
  for (std::size_t row = 0; row != n * n; ++row) {
    if (coin(rng) >= 3) continue;
    const int k = count(rng);
    const std::size_t c0 = column(rng);
    for (int j = 0; j != k; ++j) {
      const real_t v = sign(rng) ? value(rng) : -value(rng);  // never 0
      B(v, row, (c0 + j) % n);
    }
  }
  return B;
}

bool exported_operator_exists(std::size_t n, std::size_t L, std::string const &key) {
  const std::string path = "../operator/V_" + std::to_string(n * n * L * L) +
                           "_N_" + std::to_string(n) + "_L_" +
                           std::to_string(L) + "/" + key + ".ndim";
  return std::filesystem::exists(path);
}

}  // namespace

int main(int argc, char **argv) {
  const std::size_t n = argc > 1 ? std::stoul(argv[1]) : 8;
  const std::size_t L = argc > 2 ? std::stoul(argv[2]) : 3;
  if (n < 3 || L < 1) {
    std::fprintf(stderr, "need n >= 3 and L >= 1\n");
    return 2;
  }

  // ---- Problem setup (as in main.cpp) ------------------------------------
  components sbp{n, 1. / static_cast<real_t>(L)};
  sbp.n_blocks = L * L;
  sbp.n_blocks_dim = L;
  vv<std::size_t> interfaces;
  sbp.n_interfaces = make_connectivity(interfaces, L);

  auto block_grid = range_t(0, 1., L + 1);
  std::vector<range_t> grids;
  for (std::size_t i = 0; i != block_grid.size() - 1; ++i)
    grids.push_back(range_t(*block_grid[i], *block_grid[i + 1], n));

  // Refer to the static members through the type, not an instance: nvcc
  // rewrites instance.static_member into a comma expression.
  using exact = data<real_t>;
  real_t *boundary_solution = nullptr;
  compute_boundary_solution(&boundary_solution, grids,
                            boundary_functions{exact::w, exact::e, exact::s, exact::n},
                            {0., 1., 0., 1.});
  real_t *sources = nullptr;
  compute_sources(&sources, grids, exact::source_function);

  vv<std::size_t> boundary_data_map, boundary_order_map;
  make_boundary_maps(boundary_data_map, boundary_order_map, L);

  // ---- Boundary operators ------------------------------------------------
  std::vector<csr_t> B(8);
  const bool exported = exported_operator_exists(n, L, "B0");
  if (exported) {
    for (std::size_t i = 0; i != 8; ++i)
      load_operator(B[i], "B" + std::to_string(i), n, L);
  } else {
    std::mt19937_64 rng(1234);
    for (auto &b : B) b = random_boundary_operator(n, rng);
  }
  std::printf("compute_g test: n = %zu, L = %zu (%zu blocks), B operators: %s\n",
              n, L, sbp.n_blocks,
              exported ? "exported from MKL" : "random (none exported yet)");

  // ---- Host reference for every block (global numbering) -----------------
  const std::size_t n2 = n * n;
  const std::size_t face_size = n * L;
  std::vector<real_t> ref(n2 * sbp.n_blocks, 0.0);
  for (std::size_t blk = 0; blk != sbp.n_blocks; ++blk) {
    real_t *gi = &ref[n2 * blk];
    for (std::size_t face = 0; face != 4; ++face) {
      const std::size_t type = boundary_order_map[blk][face];
      if (type == 0) continue;
      const std::size_t ti = face + (type - 1) * 4;
      const std::size_t di = boundary_data_map[blk][face] - 1;
      spmv_host(B[ti], boundary_solution + face * face_size + di * n, gi, 1.0);
    }
    spmv_host(sbp.hl, sources + n2 * blk, gi, -1.0);
  }
  real_t ref_max = 0;
  for (real_t x : ref) ref_max = std::max(ref_max, std::fabs(x));

  // ---- Device inputs (shared by every simulated rank) --------------------
  real_t *d_solutions = nullptr, *d_sources = nullptr;
  CUDA_CHECK(cudaMalloc(&d_solutions, sizeof(real_t) * 4 * face_size));
  CUDA_CHECK(cudaMalloc(&d_sources, sizeof(real_t) * n2 * sbp.n_blocks));
  CUDA_CHECK(cudaMemcpy(d_solutions, boundary_solution,
                        sizeof(real_t) * 4 * face_size, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_sources, sources, sizeof(real_t) * n2 * sbp.n_blocks,
                        cudaMemcpyHostToDevice));

  std::vector<device_csr> d_B;
  for (auto const &b : B) d_B.push_back(to_device(b));
  device_csr d_h = to_device(sbp.hl);

  // ---- Run and compare ---------------------------------------------------
  constexpr real_t tolerance = 1e-13;
  bool all_pass = true;
  std::vector<real_t> g;
  for (std::size_t n_ranks = 1; n_ranks <= 3; ++n_ranks) {
    for (std::size_t rank = 0; rank != n_ranks; ++rank) {
      assign_blocks(sbp, n_ranks, rank);
      if (sbp.rank_limit_u == 0) continue;

      real_t *d_g = nullptr;
      CUDA_CHECK(cudaMalloc(&d_g, sizeof(real_t) * n2 * sbp.rank_limit_u));
      compute_g(d_g, d_B, d_h, d_solutions, d_sources, boundary_order_map,
                boundary_data_map, sbp);
      g.assign(n2 * sbp.rank_limit_u, 0.0);
      CUDA_CHECK(cudaMemcpy(g.data(), d_g, sizeof(real_t) * g.size(),
                            cudaMemcpyDeviceToHost));  // also syncs
      CUDA_CHECK(cudaFree(d_g));

      real_t err = 0;
      for (std::size_t b = 0; b != sbp.rank_limit_u; ++b) {
        const std::size_t global = sbp.rank_index_u[b];
        for (std::size_t i = 0; i != n2; ++i)
          err = std::max(err, std::fabs(g[n2 * b + i] - ref[n2 * global + i]));
      }
      err /= (ref_max > 0 ? ref_max : 1);
      const bool pass = err <= tolerance;
      all_pass = all_pass && pass;
      std::printf("  ranks %zu, rank %zu: blocks %zu..%zu  max rel err %.2e  %s\n",
                  n_ranks, rank, sbp.rank_index_u.front(),
                  sbp.rank_index_u.back(), err, pass ? "PASS" : "FAIL");
    }
  }

  for (auto &d : d_B) free_device(d);
  free_device(d_h);
  CUDA_CHECK(cudaFree(d_solutions));
  CUDA_CHECK(cudaFree(d_sources));
  free(boundary_solution);
  free(sources);

  std::printf("%s\n", all_pass ? "PASS" : "FAIL");
  return all_pass ? 0 : 1;
}
