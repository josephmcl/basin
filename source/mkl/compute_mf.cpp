#include "compute_mf.h"

// Given a vector of matrices M, and a 2-D vector of vectors F, solve 
// Mi x = Fi for all Mi in M and Fi in F. Store the result in X, aligned
// M-index major and F-index minor. 
void compute_mf(
  std::vector<real_t *>        &x,
  std::vector<sparse_matrix_t> &m,
  std::vector<real_t *>        &f,
  components             const &sbp) {
  
  matrix_descr da;
  da.type = SPARSE_MATRIX_TYPE_GENERAL;
  da.mode = SPARSE_FILL_MODE_UPPER;
	da.diag = SPARSE_DIAG_NON_UNIT;
  
  std::size_t block_index, factor_index;

  std::size_t const limit = m.size() * f.size();
  x.resize(limit);
  for (std::size_t index = 0; index != limit; ++index) {
    x[index] = (real_t *) MKL_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.n, 64);
  }

  // #pragma omp parallel for num_threads(sbp.n_threads) private(block_index, factor_index, slice_index, x_index, thread_index)
  for (std::size_t index = 0; index != limit; ++index) {
 
    block_index  = index / f.size();
  	factor_index = index % f.size();

    auto status = mkl_sparse_d_trsm (
      SPARSE_OPERATION_NON_TRANSPOSE, 1., 
      m[block_index], da, SPARSE_LAYOUT_ROW_MAJOR, 
      f[factor_index], sbp.n, CblasRowMajor, 
      x[index], CblasRowMajor);
    mkl_sparse_status(status);

    // THIS IS ONLY IMPLEMENTED FOR TRIANGULAR MATRICES USE PARADISIO

  }
}