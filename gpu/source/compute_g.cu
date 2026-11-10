
#include <cuda_runtime.h>

// make sure to cap for warp size = 32 CUDA
template <int THREADS_PER_ROW>
static __global__ void vector_csr_kernel(const int m,
                                  const int *__restrict__ row_offsets,
                                  const int *__restrict__ cols,
                                  const double *__restrict__ vals,
                                  const double *__restrict__ x,
                                  double *__restrict__ y,
                                  const double alpha,
                                  const double beta)
{
  const int row = threadIdx.y + blockDim.y * blockIdx.x;
  // CUDA *_sync shuffles need an explicit mask of participating lanes;
  // compute it before the divergent branch so every lane in the warp votes.
  const unsigned mask = __ballot_sync(0xffffffffu, row < m);
  if (row < m) {
    // determine the start and ends of each row
    int p = row_offsets[row];
    int q = row_offsets[row+1];

    // start the sparse row * vector dot product operation
    double sum = 0;
    for (int i = p + threadIdx.x; i < q; i += THREADS_PER_ROW) 
    {
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
