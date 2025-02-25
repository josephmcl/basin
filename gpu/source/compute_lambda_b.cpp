
#include "compute_lambda_b.h"

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

__global__ static void scalar_csr_kernel(const int m,
  const int *__restrict__ row_offsets,
  const int *__restrict__ cols,
  const double *__restrict__ vals,
  const double *__restrict__ x,
  double *__restrict__ y,
  const double alpha,
  const double beta)
{
  const int row = threadIdx.x + blockDim.x * blockIdx.x;
  if (row < m) {
    // determine the start and ends of each row
    int p = row_offsets[row];
    int q = row_offsets[row+1];
    // run the full sparse row * vector dot product operation
    double sum = 0;
    for (int i = p; i < q; i++) {
      sum += vals[i] * x[cols[i]];
    }
    // write to memory
    if (beta == 0) {
      y[row] = alpha * sum;
    } else {
      y[row] = alpha * sum + beta * y[row];
    }
  }
}

static void scalar_csr(
  int m,
  int threads_per_block,
  int * row_offsets,
  int * cols,
  double * vals,
  double * x,
  double * y,
  double alpha,
  double beta)
{
  int num_blocks = (m + threads_per_block - 1) / threads_per_block;
  dim3 grid(num_blocks, 1, 1);
  dim3 block(threads_per_block, 1, 1);
  scalar_csr_kernel<<<grid, block>>>(m, row_offsets, cols, vals, x, y, alpha, beta);
}

void compute_lambda_b(
    real_t             *LAMBDAb, 
    std::vector<csr_t> &Fsparse, 
    real_t             *Mg, 
    vv<std::size_t>    &FT_symbols, 
    components         &sbp) {


    std::size_t findex, j_global;
    double * lb, *mg;
    //#pragma omp parallel private(findex, lb, mg) 
    //# pragma omp for collapse(2) 
    for (std::size_t i = 0; i < sbp.n_interfaces; ++i) {
        for (std::size_t j = 0; j < sbp.rank_limit_u; ++j) {
            j_global = sbp.rank_index_u[j];
            if (FT_symbols[i][j_global] > 0) {
                findex = FT_symbols[i][j_global] - 1;
                lb = &LAMBDAb[i * sbp.n];
                // NOTE: j is the local (rank) index. j_global is the 
                //       global element index. we use j below because 
                //       Mg is a rank computation. But we j_global 
                //       elsewhere to look up where the global index 
                //       in other the symbolic data structures.
                mg = &Mg[j * sbp.n * sbp.n];
                
                vector_csr(
                    Fsparse[findex].n,
                    256,
                    64, 
                    Fsparse[findex].row_index_data(),
                    Fsparse[findex].col_index_data(),
                    Fsparse[findex].nnz(),
                    Fsparse[findex].val_data(), 
                    mg, 
                    lb, -1., 1.);
                
               /*
                scalar_csr(
                  Fsparse[findex].n,
                  // 256,
                  64, 
                  Fsparse[findex].row_index_data(),
                  Fsparse[findex].col_index_data(),
                  // Fsparse[findex].nnz(),
                  Fsparse[findex].val_data(), 
                  mg, 
                  lb, -1., 1.);
                */
                hipDeviceSynchronize();
            }
        }
    }
}


