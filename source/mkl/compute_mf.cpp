#include "compute_mf.h"

// Given a vector of matrices M, and a 2-D vector of vectors F, solve 
// Mi x = Fi for all Mi in M and Fi in F. Store the result in X, aligned
// M-index major and F-index minor. 
void compute_mf(
  std::vector<real_t *>        &x,
  vv<sparse_matrix_t> &m,
  std::vector<real_t *>        &f,
  components             const &sbp) {
  
  sparse_status_t status;

  std::size_t l;
  #pragma omp parallel for private(l) collapse(2) num_threads(sbp.n_threads)
  for (std::size_t i = 0; i != m[0].size(); ++i) {
    for (std::size_t j = 0; j != f.size(); ++j) {
      for (std::size_t k = 0; k != sbp.n; ++k) {
        auto td = omp_get_thread_num();
        l = i * f.size() + j;
        status = mkl_sparse_d_qr_solve(
          SPARSE_OPERATION_NON_TRANSPOSE, m[td][i], nullptr,
          SPARSE_LAYOUT_COLUMN_MAJOR, 1, &x[l][sbp.n * sbp.n * k] , sbp.n, 
          &f[j][sbp.n * sbp.n * k], sbp.n);
        mkl_sparse_status(status);
      }
    }
  }
}

/*
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

  // todo: give up on pardiso for now and try the QR solver

  MKL_INT iparm[64];
  iparm[0] = 1;          No solver default 
  iparm[1] = 0;          Fill-in reordering from METIS 
  iparm[3] = 0;          No iterative-direct algorithm 
  iparm[4] = 0;          No user fill-in reducing permutation 
  iparm[5] = 0;          Write solution into x 
  iparm[6] = 0;          Not in use 
  iparm[7] = 0;          Max numbers of iterative refinement steps 
  iparm[8] = 0;          Not in use 
  iparm[9] = 13;         Perturb the pivot elements with 1E-13 
  iparm[10] = 1;         Use nonsymmetric permutation and scaling MPS 
  iparm[11] = 0;         Not in use 
  iparm[12] = 0;         Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy 
  iparm[13] = 0;         Output: Number of perturbed pivots 
  iparm[14] = 0;         Not in use 
  iparm[15] = 0;         Not in use 
  iparm[16] = 0;         Not in use 
  iparm[17] = -1;        Output: Number of nonzeros in the factor LU 
  iparm[18] = -1;        Output: Mflops for LU factorization 
  iparm[19] = 0;         Output: Numbers of CG Iterations 
  
  iparm[34] = 1; 
  iparm[26] = 1;

  MKL_INT max_factors = 1;
  MKL_INT factor_choice = 1;
  MKL_INT matrix_type = 11;
  MKL_INT phase;
  MKL_INT n = static_cast<MKL_INT>(sbp.n * sbp.n);
  MKL_INT nrhs = static_cast<MKL_INT>(sbp.n);

  MKL_INT iunused;
  real_t dunused;
  MKL_INT msglvl = 1;
  MKL_INT error = 0;


  // #pragma omp parallel for num_threads(sbp.n_threads) private(block_index, factor_index, slice_index, x_index, thread_index)
  for (std::size_t index = 0; index != limit; ++index) {

    void *pt[64];
    for (int i = 0; i < 64; i++ ) {
      pt[i] = 0;
    }
 
    block_index  = index / f.size();
  	factor_index = index % f.size();

    //auto status = mkl_sparse_d_trsm (
    //  SPARSE_OPERATION_NON_TRANSPOSE, 1., 
    //  m[block_index], da, SPARSE_LAYOUT_ROW_MAJOR, 
    //  f[factor_index], sbp.n, CblasRowMajor, 
    //  x[index], CblasRowMajor);
    //mkl_sparse_status(status);

    // sparse_status_t mkl_sparse_d_export_csr (const sparse_matrix_t source, sparse_index_base_t *indexing, MKL_INT *rows, MKL_INT *cols, MKL_INT **rows_start, MKL_INT **rows_end, MKL_INT **col_indx, double **values);

    sparse_index_base_t indexing;
    MKL_INT rows, cols, *rowst, *rowe, *coli, *ia;
    real_t *vals;
    auto status = mkl_sparse_d_export_csr(
      m[block_index], &indexing, &rows, &cols, &rowst, &rowe, &coli, &vals);
    mkl_sparse_status(status);


    ia = (MKL_INT *) MKL_malloc(sizeof(MKL_INT) * ((sbp.n * sbp.n) + 1), 64);
    std::memcpy(&ia[0], &rowst[0], sizeof(MKL_INT) * sbp.n * sbp.n);
    ia[sbp.n * sbp.n] = rowe[sbp.n * sbp.n - 1];

    for (int i = 0; i < 96; ++i) {
      std::cout << vals[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < 96; ++i) {
      std::cout << coli[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < sbp.n * sbp.n + 1; ++i) {
      std::cout << ia[i] << " ";
    }
    std::cout << std::endl;

    phase = 11;
    pardiso(pt, 
      &max_factors,     Maximum # of allocd factors = 1 
      &factor_choice,   Which factor to choose = 1      
      &matrix_type,     Matrix type nonsym, posdef = 11 
      &phase, 
      &n,               Rows in A = sbp.n * sbp.n       
      vals, ia, coli,
      &iunused, 
      &nrhs,            Columns in F = sbp.n            
      iparm,
      &msglvl,
      &dunused,
      &dunused,
      &error);
    if (error != 0) {
      std::cout << "11: " << error << std::endl;
    }

    std::cout << iparm[14] << std::endl;
    std::cout << iparm[15] << std::endl;
    std::cout << iparm[16] << std::endl;
    std::cout << iparm[17] << std::endl;
    

    phase = 22;
    pardiso(pt, 
      &max_factors,     Maximum # of allocd factors = 1 
      &factor_choice,   Which factor to choose = 1      
      &matrix_type,     Matrix type nonsym, posdef = 11 
      &phase, 
      &n,               Rows in A = sbp.n * sbp.n       
      vals, ia, coli,
      &iunused, 
      &nrhs,            Columns in F = sbp.n            
      iparm,
      &msglvl,
      &dunused,
      &dunused,
      &error);
    if (error != 0) {
      std::cout << "22: " << error << std::endl;
    }

    x[0] = (real_t *) MKL_malloc(sizeof(real_t) * sbp.n * sbp.n * sbp.n, 64);

    phase = 33;
    pardiso(pt, 
      &max_factors,     Maximum # of allocd factors = 1 
      &factor_choice,   Which factor to choose = 1      
      &matrix_type,     Matrix type nonsym, posdef = 11 
      &phase, 
      &n,               Rows in A = sbp.n * sbp.n       
      vals, ia, coli,
      &iunused, 
      &nrhs,            Columns in F = sbp.n            
      iparm,
      &msglvl,
      (void *)f[factor_index],
      (void *)x[0],
      &error);
    if (error != 0) {
      std::cout << "33: " << error << std::endl;
    }
    MKL_free(rowst);
    MKL_free(rowe);
    MKL_free(coli);
    MKL_free(vals);
    MKL_free(ia); */
