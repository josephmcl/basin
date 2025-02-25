#include "compute_rhs.h"


#include "hip/hip_runtime.h"

static unsigned long prevPowerOf2(unsigned long v) {
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v >> 1;
}

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
  if (row < m) {
    // determine the start and ends of each row
    int p = row_offsets[row];
    int q = row_offsets[row+1];

    // start the sparse row * vector dot product operation
    double sum = 0;
    for (int i = p + threadIdx.x; i < q; i += THREADS_PER_ROW) {
      sum += vals[i] * x[cols[i]];
    }

    // finish the sparse row * vector dot product operation
#pragma unroll
    for (int i = THREADS_PER_ROW >> 1; i > 0; i >>= 1)
      sum += __shfl_down(sum, i, THREADS_PER_ROW);

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

static void vector_csr(int m,
                int threads_per_block,
                int warpSize,
                int * row_offsets,
                int * cols,
                int nnz,
                double * vals, 
                double * x,
                double * y,
                double alpha,
                double beta)
{

  int nnz_per_row = nnz / m;
  int threads_per_row = prevPowerOf2(nnz_per_row);
  // limit the number of threads per row to be no larger than the wavefront (warp) size
  threads_per_row = threads_per_row > warpSize ? warpSize : threads_per_row;


  //std::cout << m << " " << nnz << std::endl;
  //std::cout << threads_per_block << " " << threads_per_row << std::endl;

  int rows_per_block;
  if (threads_per_row > 0) {
    rows_per_block = threads_per_block / threads_per_row;
  }
  else {
    rows_per_block = 1;
  }
  int num_blocks = (m + rows_per_block - 1) / rows_per_block;

  dim3 grid(num_blocks, 1, 1);
  dim3 block(threads_per_row, rows_per_block, 1);
  
  if (threads_per_row <= 2)
      vector_csr_kernel<2><<<grid, block>>>(m, row_offsets, cols, vals, x, y, alpha, beta);
  else if (threads_per_row <= 4)
      vector_csr_kernel<4><<<grid, block>>>(m, row_offsets, cols, vals, x, y, alpha, beta);
  else if (threads_per_row <= 8)
      vector_csr_kernel<8><<<grid, block>>>(m, row_offsets, cols, vals, x, y, alpha, beta);
  else if (threads_per_row <= 16)
      vector_csr_kernel<16><<<grid, block>>>(m, row_offsets, cols, vals, x, y, alpha, beta);
  else if (threads_per_row <= 32)
      vector_csr_kernel<32><<<grid, block>>>(m, row_offsets, cols, vals, x, y, alpha, beta);
  else
      vector_csr_kernel<64><<<grid, block>>>(m, row_offsets, cols, vals, x, y, alpha, beta);
    
}

void compute_rhs(
    real_t *rhs,
    real_t *g, 
    std::vector<csr_t> &F,
    real_t *lambda,
    vv<std::size_t> &F_symbols,
    components &sbp) {

    // Compute rhs = g - F * lambda

    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> temp;
    for (auto &i : sbp.rank_index_u) {
        for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
            if (F_symbols[i][j] != 0) {
                temp.push_back({i - sbp.rank_index_u[0], j, F_symbols[i][j] - 1});
            }
        }
    }


    std::size_t i, j, k;
    double *l, *r;
    // #pragma omp parallel for private(i, j, k, l, r) 
    for (std::size_t a = 0; a < temp.size(); ++a) {
        i = std::get<0>(temp[a]);
        j = std::get<1>(temp[a]);
        k = std::get<2>(temp[a]);
        l = &lambda[j * sbp.n];
        r = &rhs[i * sbp.n * sbp.n];
        
        /*
        status = mkl_sparse_d_mv(
            SPARSE_OPERATION_TRANSPOSE, -1., 
            F[k], da,
            l, 
            1., r);
        mkl_sparse_status(status); 
        */

        vector_csr(
            F[k].n,
            256,
            64, 
            F[k].row_index_data(),
            F[k].col_index_data(),
            F[k].nnz(),
            F[k].val_data(), 
            l, 
            r, 1., 1.);
    }
    
    for (std::size_t i = 0; i != sbp.n * sbp.n * sbp.rank_limit_u; ++i) {
        rhs[i] += g[i];
    }

} 