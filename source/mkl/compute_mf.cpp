#include "compute_mf.h"

std::size_t compute_mf_pre(
  std::vector<real_t *>        &x,
  vv<sparse_matrix_t>          &m,
  std::vector<real_t *>        &f,
  components             const &sbp) {
  
  x.resize(m[0].size() * f.size());
  std::size_t sz;
  for (std::size_t index = 0; index != x.size(); ++index) {
      sz = sbp.n * sbp.n * sbp.n;
      x[index] = (real_t *) mkl_malloc(sizeof(real_t) * sz, 64);
      memset(x[index], 0, sizeof(real_t) * sz);
  }
  return 0;
}


// Given a vector of matrices M, and a 2-D vector of vectors F, solve 
// Mi x = Fi for all Mi in M and Fi in F. Store the result in X, aligned
// M-index major and F-index minor. 
std::size_t compute_mf_call(
  std::vector<real_t *>        &x,
  vv<sparse_matrix_t>          &m,
  std::vector<real_t *>        &f,
  components             const &sbp) {
  
  sparse_status_t status;

  std::size_t l;
  #pragma omp parallel for private(l) collapse(2) 
  for (std::size_t i = 0; i < m[0].size(); ++i) {
    for (std::size_t j = 0; j < f.size(); ++j) {
      for (std::size_t k = 0; k < sbp.n; ++k) {
        auto td = omp_get_thread_num();
        l = i * f.size() + j;
        status = mkl_sparse_d_qr_solve(
          SPARSE_OPERATION_NON_TRANSPOSE, 
          m[td][i], 
          nullptr,
          SPARSE_LAYOUT_COLUMN_MAJOR, 
          1, 
          &x[l][sbp.n * sbp.n * k] , 
          sbp.n, 
          &f[j][sbp.n * sbp.n * k], 
          sbp.n);
        mkl_sparse_status(status);
      }
    }
  }
  return 0;
}

mf_instrument compute_mf(compute_mf_pre, compute_mf_call, std::nullopt);
