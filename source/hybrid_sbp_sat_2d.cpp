#include "hybrid_sbp_sat_2d.h"


void sbp_sat::x2::
petsc_hybridized_poisson(sbp_sat::real_v             &result,
                         sbp_sat::domain_v     const &domain, 
                         sbp_sat::block_v      const &blocks, 
                         sbp_sat::boundary_vx2 const &boundaries) {

   
  using namespace sbp_sat;

  nat_t local_problems = blocks.size() * blocks.size();

  /*  Create petsc objects for the block diagonal matrix, M. 

          [ m_0           ] 
      M = [     ...       ] in [ M  F ]  
          [         m_n-1 ]    [ Ft D ]   

      From the domain blocks 

      | b0  | ... | ...  |
      |-----|-----|------|
      | ... | ... | ...  |
      |-----|-----|------|
      | ... | ... | bn-1 |  */


  auto M = std::vector<petsc_matrix>(local_problems); 

  std::cout << local_problems << " local problems." << std::endl;

  for (std::size_t row = 0; row != blocks.size(); ++row) {
    for (std::size_t col = 0; col != blocks.size(); ++col) {

      auto matrix_index = blocks.size() * row + col;
      auto matrix_block = M[matrix_index];

      auto [row_x1, col_x1] = blocks[row];
      auto [row_x2, col_x2] = blocks[col]; 

      (void) col_x1; (void) row_x2; /* unused */

      
      MatCreateSeqAIJ(PETSC_COMM_SELF, row_x1.size(), col_x2.size(), 
        4, nullptr, &matrix_block);

      x2::write_m(matrix_block, blocks[row], blocks[col], 
        boundaries[row], boundaries[col]);

      MatAssemblyBegin(matrix_block, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(matrix_block, MAT_FINAL_ASSEMBLY);

    }
  }

  /* Free all petsc objects. */
  for (auto &e : M) MatDestroy(&e);

}

void sbp_sat::x2::
write_m(petsc_matrix       &m,
        block_t      const &block_x1, 
        block_t      const &block_x2, 
        boundary_tx2 const &boundary_x1, 
        boundary_tx2 const &boundary_x2) {

  auto [row_x1, col_x1] = block_x1;
  auto [row_x2, col_x2] = block_x2;

  (void) col_x1; (void) row_x2; // unused;

  auto x1_spacing = row_x1.to - row_x1.from;
  auto x2_spacing = col_x2.to - col_x2.from;

  auto x1_sz = row_x1.size() - 1;
  auto x2_sz = col_x2.size() - 1;

  auto x1_spacing_square = x1_spacing * x1_spacing;
  auto x2_spacing_square = x2_spacing * x2_spacing;


  auto hx1 = numerical::operators::H(row_x1.size(), 2, 0, x1_spacing);
  auto hx2 = numerical::operators::H(col_x2.size(), 2, 0, x2_spacing);

  std::cout << x1_spacing << " " << x2_spacing << std::endl;
  std::cout << x1_spacing_square << " " << x2_spacing_square << std::endl;

  petsc_matrix d2x1;
  petsc_matrix d2x2;

  MatCreateSeqAIJ(PETSC_COMM_SELF, row_x1.size(), row_x1.size(), 4, 
    nullptr, &d2x1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, col_x2.size(), col_x2.size(), 4, 
    nullptr, &d2x2);

  write_d2_h1(d2x1, hx1, row_x1, x1_spacing_square / (x1_sz * x1_sz));
  write_d2_h1(d2x2, hx2, col_x2, x2_spacing_square / (x2_sz * x2_sz));

  MatAssemblyBegin(d2x1, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(d2x1, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(d2x2, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(d2x2, MAT_FINAL_ASSEMBLY);

  MatView(d2x1, PETSC_VIEWER_STDOUT_WORLD);
  MatView(d2x2, PETSC_VIEWER_STDOUT_WORLD);

  MatDestroy(&d2x1);
  MatDestroy(&d2x2);

  // MatSeqAIJKron(x1, x2, MAT_REUSE_MATRIX , &m)

}

void sbp_sat::x2::write_d2_h1(
  petsc_matrix                   &M, 
  std::vector<long double> const &h, 
  range_t                  const &local,
  long double              const  spacing_square) {
    
  /* Initialize the first local skew row. */
  MatSetValue(M, 0, 0,  1. / spacing_square * h[0], ADD_VALUES);
  MatSetValue(M, 0, 1, -2. / spacing_square * h[0], ADD_VALUES);
  MatSetValue(M, 0, 2,  1. / spacing_square * h[0], ADD_VALUES); 

  /* Initialize the final local skew row. */
  auto n = local.size() - 1;
  MatSetValue(M, n, n - 2,  1. / spacing_square * h[n], ADD_VALUES);
  MatSetValue(M, n, n - 1, -2. / spacing_square * h[n], ADD_VALUES);
  MatSetValue(M, n, n,      1. / spacing_square * h[n], ADD_VALUES); 

  /* Initialize the interior local diagonal rows. */
  for (auto it = local.begin() + 1; it != local.end() - 1; ++it) {  
    auto i = it.index;
    MatSetValue(M, i, i - 1,  1. / spacing_square * h[i], ADD_VALUES);
    MatSetValue(M, i, i,     -2. / spacing_square * h[i], ADD_VALUES);
    MatSetValue(M, i, i + 1,  1. / spacing_square * h[i], ADD_VALUES);  
  }
}