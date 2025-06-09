
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
    double *lb, *mg;
    //#pragma omp parallel private(findex, lb, mg) 
    //# pragma omp for collapse(2) 
    //#pragma omp parallel for collapse(2) private(findex, lb, mg) 
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
               
                /*
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
                */
                
                
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
                hipDeviceSynchronize();
            }
        }
    }
}

void compute_lambda_b_csr(
  real_t     *λb, 
  csr_t      *F, 
  real_t     *Mg, 
  components &sbp) {

  /*
  vector_csr(
    F->n,
    256,
    64, 
    F->row_index_data(),
    F->col_index_data(),
    F->nnz(),
    F->val_data(), 
    Mg, 
    λb, -1., 1.);*/
  
  
  scalar_csr(
    F->n,
    // 256,
    64, 
    F->row_index_data(),
    F->col_index_data(),
    // Fsparse[findex].nnz(),
    F->val_data(), 
    Mg, 
    λb, -1., 1.);

  hipDeviceSynchronize();

}




// 1.359493 1.397719 0.472190 0.120601 1.370889 1.184983 0.209268 0.043572 0.079335 0.265526 0.210132 -0.013702 0.013702 -0.210132 -0.265526 -0.079335 -0.043572 -0.209268 -1.184983 -1.370889 -0.120601 -0.472190 -1.397719 -1.359493 0.041010 0.412427 0.437594 0.072205 0.043500 0.339957 0.287821 -0.046427 -0.005432 -0.230141 -0.568669 -0.484457 0.484457 0.568669 0.230141 0.005432 0.046427 -0.287821 -0.339957 -0.043500 -0.072205 -0.437594 -0.412427 -0.041010 
// 2.878014 1.765460 -0.245490 0.330811 -0.563814 -0.810516 1.906716 2.878014 0.000000 -0.141257 0.565026 0.894625 -0.894625 -0.565026 0.141257 0.000000 -2.878014 -1.906716 0.810516 0.563814 -0.330811 0.245490 -1.765460 -2.878014 -0.016341 1.113706 1.170839 0.110032 1.362638 1.221968 0.246907 0.036126 -0.025716 -1.394232 -1.425698 -0.091704 0.091704 1.425698 1.394232 0.025716 -0.036126 -0.246907 -1.221968 -1.362638 -0.110032 -1.170839 -1.113706 0.016341 

//2.878014 1.765460 -0.245490 0.330811 -0.563814 -0.810516 1.906716 2.878014 0.000000 -0.141257 0.565026 0.894625 -0.894625 -0.565026 0.141257 0.000000 -2.878014 -1.906716 0.810516 0.563814 -0.330811 0.245490 -1.765460 -2.878014 -0.016341 1.113706 1.170839 0.110032 1.362638 1.221968 0.246907 0.036126 -0.025716 -1.394232 -1.425698 -0.091704 0.091704 1.425698 1.394232 0.025716 -0.036126 -0.246907 -1.221968 -1.362638 -0.110032 -1.170839 -1.113706 0.016341