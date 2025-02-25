
#include "compute_g.h"

#include "hip/hip_runtime.h"

#define HIP_ENABLE_PRINTF

#define HIP_CHECK(command)                                    \
{                                                             \
  hipError_t stat = (command);                                \
  if(stat != hipSuccess)                                      \
  {                                                           \
    std::cerr << "HIP error: " << hipGetErrorString(stat) <<  \
    " in file " << __FILE__ << ":" << __LINE__ << std::endl;  \
    exit(-1);                                                 \
  }                                                           \
}

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
  std::cout << nnz << " " << m << std::endl;
  int threads_per_row = prevPowerOf2(nnz_per_row);
  // limit the number of threads per row to be no larger than the wavefront (warp) size
  threads_per_row = threads_per_row > warpSize ? warpSize : threads_per_row;
  int rows_per_block = threads_per_block / threads_per_row;
  int num_blocks = (m + rows_per_block - 1) / rows_per_block;

  dim3 grid(num_blocks, 1, 1);
  dim3 block(threads_per_row, rows_per_block, 1);

  std::cout << prevPowerOf2(nnz_per_row) << " " << std::endl;
  std::cout << "dims " << num_blocks << " " << threads_per_row << " " << rows_per_block << std::endl;
  
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

void compute_g(
    real_t                  **g, 
    std::vector<csr_t>       &boundaries,
    real_t                   *solutions,
    real_t                   *sources, 
    vv<std::size_t>          &boundary_type_map,
    vv<std::size_t>          &boundary_data_map,
    components               &sbp) {


    // gi = [n2 vec], solution = [n vec], B[i] = [n2 x n] mat


    // compute ith: (boundary * solution ...) - H_tilde sources[i]
    //                       bm       bm         stat     m2 block 
    //    v         v        v        v         v         v      
    // (b_LB_W * g_LB_W + b_LB_S * g_LB_S) - H_tilde * F_LB[:]


    (*g) = (real_t *) malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u);
    memset(*g, 0, sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u);

    real_t *gtemp = (real_t *) malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u);
    memset(gtemp, 0, sizeof(real_t) * sbp.n * sbp.n * sbp.rank_limit_u);

    std::size_t n2 = sbp.n * sbp.n;
    std::size_t face_size = sbp.n * sbp.n_blocks_dim;
    real_t *gi, *gti, *solution, *source;


    constexpr double alpha = 1.;

    std::size_t relblock; 
    for (std::size_t block = 0; block != sbp.rank_limit_u; ++block) {
        
        
      //std::cout << "block " << block << std::endl;
      relblock = sbp.rank_index_u[block];
      // block is the local index local -> 0, 1, 2, ... 
      // but relblock is element index --> k, k+1, k+2, ... 
      gi = &(*g)[0] + (n2 * relblock); 
      gti = &(gtemp)[0] + (n2 * relblock);
      
      // 4 is known quantity as our blocks are rect. and orth.
      constexpr std::size_t faces = 4;
      for (std::size_t face = 0; face < faces; ++face) {

          
        // b_type_index gives the types (direction) of a boundary
        auto b_type_index = boundary_type_map[relblock][face];
        if (b_type_index != 0) {
        
          // std::cout << "face - " << face + ((b_type_index - 1) * faces) << std::endl;
          // ti gives face + b_type_index 
          auto ti = face + ((b_type_index - 1) * faces);
          auto boundary = boundaries[ti];

          // di gives which index of boundary data on a given face 
          // [ 0][ 1][ 2] 
          // each index is length n
          // [..face 0..][..face 1..][..face 2..][..face 3..] 
          // each face is length n * n_block_dim
          // [...............boundary solution..............]
          auto di = boundary_data_map[relblock][face] - 1;
          solution = &solutions[0] + (face * face_size) + (di * sbp.n);
          // gi += boundary:matrix * solution:vector
          
          scalar_csr(
              boundary.n,
              32, // 64, 
              boundary.row_index_data(),
              boundary.col_index_data(),
              // boundary.nnz(),
              boundary.val_data(), 
              solution, 
              gi, 1., 0.);
          hipDeviceSynchronize();
        }
      }
      
      source = &sources[0] + (sbp.n * sbp.n * relblock);

      // std::size_t source_offset = relblock * sbp.n * sbp.n;

      for (std::size_t i = 0; i != sbp.n * sbp.n; ++i) {
          gti[i] = -source[i];
      }

      scalar_csr(
        sbp.hl.n,
        //256,
        64, 
        sbp.hl.row_index_data(),
        sbp.hl.col_index_data(),
        //sbp.hl.nnz(),
        sbp.hl.val_data(), 
        gti, 
        gi, 1., 1.);
      hipDeviceSynchronize();
  }
}