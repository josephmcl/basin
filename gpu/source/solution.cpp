



// Given b as a single vector and A as a vector of solvers, 
// first move b into bs, a vector of vectors,
// then compute solve(A_i, b_i) for all A_i, b_i in A, bs
// store  

/*
void sbp_sat::x2::compute_solution(
  std::vector<petsc_vector>       &x,
  vv<petsc_solver>          const &solvers, 
  petsc_vector              const &b, 
  components                const &sbp) {

  std::size_t const n2 = sbp.n * sbp.n;

  x.resize(sbp.n_blocks);

  std::vector<int> writei; writei.resize(n2);
  for (std::size_t i = 0; i != writei.size(); ++i) {
    writei[i] = i;
  }

  std::vector<petsc_vector> bs; 
  bs.resize(sbp.n_blocks);


  // Read out all b values.
  double const *val;
  VecGetArrayRead(b, &val);

  for (std::size_t block = 0; block != sbp.n_blocks; ++block) {
    VecCreateSeq(PETSC_COMM_SELF, n2, &bs[block]);
    VecCreateSeq(PETSC_COMM_SELF, n2, &x[block]);
  }

  std::vector<std::vector<double>> bs_; bs_.resize(sbp.n_blocks);
  #pragma omp parallel for
  for (std::size_t block = 0; block != sbp.n_blocks; ++block) {
    // double const *v = val + (n2 * block);
    // std::cout << v << " " << v + n2 << std::endl;
    bs_[block].resize(n2);
    bs_[block].assign(&val[n2 * block], &val[n2 * (block + 1)]);
    //VecCreateSeq(PETSC_COMM_SELF, n2, &bs[block]);
    //VecCreateSeq(PETSC_COMM_SELF, n2, &x[block]);
  }

  //#pragma omp parallel for 

  std::size_t j, thread_index;
  #pragma omp parallel for private(j, thread_index)
  for (std::size_t i = 0; i != sbp.n_blocks_dim; ++i) {
    j = sbp.n_blocks_dim * i;
    thread_index = omp_get_thread_num();
    VecSetValues(bs[j], n2, &writei[0], &bs_[j][0], INSERT_VALUES);
    finalize<fw>(bs[j]);

    KSPSolve(solvers[0][thread_index], bs[j], x[j]); 
    finalize<fw>(x[j]);
  }

  #pragma omp parallel for private(j, thread_index)
  for (std::size_t i = 0; i != sbp.n_blocks_dim; ++i) {
    j = (sbp.n_blocks_dim * i) + sbp.n_blocks_dim - 1;
    thread_index = omp_get_thread_num();
    VecSetValues(bs[j], n2, &writei[0], &bs_[j][0], INSERT_VALUES);
    finalize<fw>(bs[j]);

    KSPSolve(solvers[2][thread_index], bs[j], x[j]); 
    finalize<fw>(x[j]);
  }

  #pragma omp parallel for private(j, thread_index)
  for (std::size_t i = 0; i != sbp.n_blocks_dim * (sbp.n_blocks_dim - 2); ++i) {
    j = i + ((i / (sbp.n_blocks_dim - 2)) * 2) + 1;
    thread_index = omp_get_thread_num();
    VecSetValues(bs[j], n2, &writei[0], &bs_[j][0], INSERT_VALUES);
    finalize<fw>(bs[j]);

    KSPSolve(solvers[1][thread_index], bs[j], x[j]); 
    finalize<fw>(x[j]);
  }

  // VecView(x[solvers.size() - 1], PETSC_VIEWER_STDOUT_SELF);

}

*/