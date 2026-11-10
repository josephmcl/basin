#include "spmv.cuh"
#include "cuda_check.cuh"

namespace {

__global__ void scalar_csr_kernel(const int m,
                                  const int *__restrict__ row_offsets,
                                  const int *__restrict__ cols,
                                  const double *__restrict__ vals,
                                  const double *__restrict__ x,
                                  double *__restrict__ y,
                                  const double alpha,
                                  const double beta) {
  const int row = threadIdx.x + blockDim.x * blockIdx.x;
  if (row < m) {
    double sum = 0;
    for (int i = row_offsets[row]; i < row_offsets[row + 1]; ++i) {
      sum += vals[i] * x[cols[i]];
    }
    if (beta == 0) {
      y[row] = alpha * sum;
    } else {
      y[row] = alpha * sum + beta * y[row];
    }
  }
}

// Launch with blockDim.x == THREADS_PER_ROW and blockDim.x * blockDim.y a
// multiple of 32, so each row's lanes sit inside one warp and every warp is
// full when __ballot_sync runs.
template <int THREADS_PER_ROW>
__global__ void vector_csr_kernel(const int m,
                                  const int *__restrict__ row_offsets,
                                  const int *__restrict__ cols,
                                  const double *__restrict__ vals,
                                  const double *__restrict__ x,
                                  double *__restrict__ y,
                                  const double alpha,
                                  const double beta) {
  const int row = threadIdx.y + blockDim.y * blockIdx.x;
  // CUDA *_sync shuffles need an explicit mask of participating lanes;
  // compute it before the divergent branch so every lane in the warp votes.
  const unsigned mask = __ballot_sync(0xffffffffu, row < m);
  if (row < m) {
    // determine the start and ends of each row
    const int p = row_offsets[row];
    const int q = row_offsets[row + 1];

    // start the sparse row * vector dot product operation
    double sum = 0;
    for (int i = p + threadIdx.x; i < q; i += THREADS_PER_ROW) {
      sum += vals[i] * x[cols[i]];
    }

    // finish the sparse row * vector dot product operation
#pragma unroll
    for (int i = THREADS_PER_ROW >> 1; i > 0; i >>= 1)
      sum += __shfl_down_sync(mask, sum, i, THREADS_PER_ROW);

    // write to memory
    if (!threadIdx.x) {
      if (beta == 0) {
        y[row] = alpha * sum;
      } else {
        y[row] = alpha * sum + beta * y[row];
      }
    }
  }
}

}  // namespace

void scalar_csr(device_csr const &A, real_t const *x, real_t *y,
                real_t alpha, real_t beta) {
  if (A.n == 0) return;
  constexpr int threads_per_block = 128;
  const int num_blocks = (A.n + threads_per_block - 1) / threads_per_block;
  scalar_csr_kernel<<<num_blocks, threads_per_block>>>(
      A.n, A.row, A.col, A.val, x, y, alpha, beta);
  CUDA_CHECK(cudaGetLastError());
}

void vector_csr(device_csr const &A, real_t const *x, real_t *y,
                real_t alpha, real_t beta) {
  if (A.n == 0) return;
  constexpr int threads_per_block = 256;

  // Largest power of two <= average nonzeros per row, clamped to [2, 32].
  // (The HIP version could produce 0 here, i.e. a launch with blockDim.x == 0,
  // whenever a matrix averaged less than one nonzero per row.)
  const int nnz_per_row = A.nnz / A.n;
  int threads_per_row = 2;
  while (threads_per_row * 2 <= nnz_per_row && threads_per_row < 32)
    threads_per_row *= 2;

  const int rows_per_block = threads_per_block / threads_per_row;
  const int num_blocks = (A.n + rows_per_block - 1) / rows_per_block;
  const dim3 grid(num_blocks);
  const dim3 block(threads_per_row, rows_per_block);

  switch (threads_per_row) {
    case 2:  vector_csr_kernel<2> <<<grid, block>>>(A.n, A.row, A.col, A.val, x, y, alpha, beta); break;
    case 4:  vector_csr_kernel<4> <<<grid, block>>>(A.n, A.row, A.col, A.val, x, y, alpha, beta); break;
    case 8:  vector_csr_kernel<8> <<<grid, block>>>(A.n, A.row, A.col, A.val, x, y, alpha, beta); break;
    case 16: vector_csr_kernel<16><<<grid, block>>>(A.n, A.row, A.col, A.val, x, y, alpha, beta); break;
    default: vector_csr_kernel<32><<<grid, block>>>(A.n, A.row, A.col, A.val, x, y, alpha, beta); break;
  }
  CUDA_CHECK(cudaGetLastError());
}
