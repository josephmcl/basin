#include "hybrid_sbp_sat_2d.h"

void sbp_sat::x2::copy_approx(
  petsc_vector &approx, 
  std::vector<petsc_vector> const &x, 
  components const &sbp) {

  // Allocate contiguous approx matrix.
  VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &approx);

  std::vector<int> index(sbp.n);
  const double *values;

  for (std::size_t block = 0; block != sbp.n_blocks; ++block) {
    for (std::size_t i = 0; i != sbp.n; ++i) {
      index[i] = block * sbp.n + i;
    }

    // Copy values from blocked vector to contiguous.
    VecGetArrayRead(x[block], &values);
    VecSetValues(approx, sbp.n, &index[0], values, INSERT_VALUES);
  }
}

void sbp_sat::x2::compute_exact(
  petsc_vector &exact, 
  components const &sbp) {

  VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &exact);

  std::vector<double> value(sbp.n * sbp.n);

  auto span = 1 / (sbp.n_blocks_dim * sbp.n - 1);

  for (std::size_t bx = 0; bx != sbp.n_blocks_dim; ++bx) {
    for (std::size_t by = 0; by != sbp.n_blocks_dim; ++by) {
      for (std::size_t nx = 0; nx != sbp.n; ++nx) {
        for (std::size_t ny = 0; ny != sbp.n; ++ny) {
            
          //analytical_solution<double>(x, y);

        }
      }
    }
  }
    
}

void sbp_sat::x2::
print_convergence(
  petsc_vector &approx, 
  petsc_vector &exact, 
  components const &sbp) {

  petsc_vector left, right, temp;
  double result = 0.;

  VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &left);
  VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &right);
  VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp);
  finalize<fw>(temp);


  VecWAXPY(left, -1., exact, approx);
  VecWAXPY(right, 1., exact, approx);

  MatMult(sbp.hl, left, temp);
  //VecTDot(temp, right, &result);

  std::cout << "... " << result << std::endl;

  // destroy<fw>(left);
  // destroy<fw>(right);
  // destroy<fw>(temp);
}