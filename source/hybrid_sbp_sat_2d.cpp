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

      break;

      MatAssemblyBegin(matrix_block, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(matrix_block, MAT_FINAL_ASSEMBLY);

    }

    break;
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

  // hard code for now 
  real_t τ = 1.;
  real_t β = 1.;

  // hard code the second order bs matrix for now. 
  vector_t bsx1 = {
    (3./2.) / x1_spacing * (x1_sz), 
        -2. / x1_spacing * (x1_sz),   // NOTE: maybe flip sz's? 
    (1./2.) / x1_spacing * (x1_sz)};

  vector_t bsx2 = {
    (3./2.) / x2_spacing * (x2_sz), 
        -2. / x2_spacing * (x2_sz), 
    (1./2.) / x2_spacing * (x2_sz)};

  matrix<fw> d2x1, d2x2, h1x1, h1x2, Ax1, Ax2;
  MatCreateSeqAIJ(PETSC_COMM_SELF, row_x1.size(), row_x1.size(), 4, 
    nullptr, &d2x1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, col_x2.size(), col_x2.size(), 4, 
    nullptr, &d2x2);
  MatCreateSeqAIJ(PETSC_COMM_SELF, row_x1.size(), row_x1.size(), 1, 
    nullptr, &h1x1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, col_x2.size(), col_x2.size(), 1, 
    nullptr, &h1x2);

  // auto A_size = row_x1.size() * col_x2.size();

  // MatCreateSeqAIJ(PETSC_COMM_SELF, A_size, A_size, 1, nullptr, &Ax1);

  write_d2_h1(d2x1, hx1, row_x1, x1_spacing_square / (x1_sz * x1_sz), -1.);
  write_d2_h1(d2x2, hx2, col_x2, x2_spacing_square / (x2_sz * x2_sz), -1.);
  write_h1(h1x1, hx1);
  write_h1(h1x2, hx2);
  finalize<fw>(d2x1);
  finalize<fw>(d2x2);
  finalize<fw>(h1x1);
  finalize<fw>(h1x2);

  MatSeqAIJKron(h1x2, d2x1, MAT_INITIAL_MATRIX, &Ax1);
  MatSeqAIJKron(d2x2, h1x1, MAT_INITIAL_MATRIX, &Ax2);

  // finalize<fw>(Ax1);
  finalize<fw>(Ax2);
  finalize<fw>(Ax1);
  
  auto sz = row_x1.size() * col_x2.size();
  MatCreate(PETSC_COMM_SELF, &m);
  MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, sz, sz);
  MatSetType(m, MATCOMPOSITE);

  MatCompositeAddMat(m, Ax1);
  MatCompositeAddMat(m, Ax2);

  // dirichlet neumann
  add_boundary<x, left,  neumann>(m, row_x1, col_x2, bsx1, hx2, β, τ);
  add_boundary<x, right, neumann>(m, row_x1, col_x2, bsx1, hx2, β, τ);
  add_boundary<y, left,  neumann>(m, row_x1, col_x2, bsx1, hx2, β, τ);
  add_boundary<y, right, neumann>(m, row_x1, col_x2, bsx1, hx2, β, τ);
  
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m,   MAT_FINAL_ASSEMBLY);
  MatCompositeSetMatStructure(m, DIFFERENT_NONZERO_PATTERN);
  MatCompositeMerge(m);
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);

  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  MatView(m, PETSC_VIEWER_STDOUT_SELF);

  MatDestroy(&d2x1);
  MatDestroy(&d2x2);
  
  MatDestroy(&h1x1);
  MatDestroy(&h1x2);

  MatDestroy(&Ax1);
  MatDestroy(&Ax2);

}

void sbp_sat::x2::write_d2_h1(
  petsc_matrix       &M, 
  real_v       const &h, 
  range_t      const &local,
  real_t       const  spacing_square, 
  real_t       const  coeff) {
    
  /* Initialize the first local skew row. */
  MatSetValue(M, 0, 0,  1. / spacing_square * h[0] * coeff, ADD_VALUES);
  MatSetValue(M, 0, 1, -2. / spacing_square * h[0] * coeff, ADD_VALUES);
  MatSetValue(M, 0, 2,  1. / spacing_square * h[0] * coeff, ADD_VALUES); 

  /* Initialize the final local skew row. */
  auto n = local.size() - 1;
  MatSetValue(M, n, n - 2,  1. / spacing_square * h[n] * coeff, ADD_VALUES);
  MatSetValue(M, n, n - 1, -2. / spacing_square * h[n] * coeff, ADD_VALUES);
  MatSetValue(M, n, n,      1. / spacing_square * h[n] * coeff, ADD_VALUES); 

  /* Initialize the interior local diagonal rows. */
  for (auto it = local.begin() + 1; it != local.end() - 1; ++it) {  
    auto i = it.index;
    MatSetValue(M, i, i - 1,  1. / spacing_square * h[i] * coeff, ADD_VALUES);
    MatSetValue(M, i, i,     -2. / spacing_square * h[i] * coeff, ADD_VALUES);
    MatSetValue(M, i, i + 1,  1. / spacing_square * h[i] * coeff, ADD_VALUES);  
  }
}

void sbp_sat::x2::write_h1(
  petsc_matrix                   &M, 
  std::vector<long double> const &h) {

  for (std::size_t i = 0; i < h.size(); ++i) {
    MatSetValue(M, i, i,  h[i], ADD_VALUES);
  }
}



void sbp_sat::x2::write_boundary_x1_left(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bsx1,
  vector_t     const &hx2, 
  real_t       const  β,
  real_t       const  τ,
  int          const  order) {

  matrix<fw> result;
  nat_t size = x1.size() * x2.size();

  if (order == 1) {
    make_local_sparse_matrix<fw>(result, size, size, 3);
    boundary_diagonal<0, 0>(result, x1.size(), x2.size(), hx2, τ);
    boundary_bs<0, 0, transposed>(result, x1.size(), x2.size(), bsx1, hx2, -β);
  }
  else { /* second order */
   
    matrix<fw> a;
    matrix<fw> b;
    make_local_sparse_matrix<fw>(a, size, size, 3);
    make_local_sparse_matrix<fw>(b, size, size, 3);
    boundary_diagonal<0, 0>(a, x1.size(), x2.size(), hx2);
    boundary_bs<0, 0, transposed>(a, x1.size(), x2.size(), bsx1, hx2, -1. / τ);
    boundary_bs<0, 0>(b, x1.size(), x2.size(), bsx1);
    finalize<fw>(a);
    finalize<fw>(b);
    matmul<fw>(a, b, result);
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    MatView(a, PETSC_VIEWER_STDOUT_SELF);
    MatView(b, PETSC_VIEWER_STDOUT_SELF);
    destroy<fw>(a);
    destroy<fw>(b);
  }
  finalize<fw>(result);
  MatCompositeAddMat(M, result);
}


void sbp_sat::x2::write_boundary_x1_right(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bsx1,
  vector_t     const &hx2, 
  real_t       const  β,
  real_t       const  τ,
  int          const  order) {

  auto size = (x1.size() * x2.size()) - 1;

  for (auto e1 = x1.begin(); e1 != x1.end(); ++e1) {
    auto i = e1.index;
    MatSetValue(M, size - i, size - i, 
      τ * hx2[i] - β * hx2[i] * bsx1[0], ADD_VALUES);
    MatSetValue(M, size - i - x1.size(), size - i,       
      -β * hx2[i] * bsx1[1], ADD_VALUES);
    MatSetValue(M, size - i - (2 * x1.size()), size - i, 
      -β * hx2[i] * bsx1[2], ADD_VALUES);
  }
}

void sbp_sat::x2::write_boundary_x2_left(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bsx2,
  vector_t     const &hx1, 
  real_t       const  β,
  real_t       const  τ,
  int          const  order) {

  auto x2_sz = x2.size();

  for (auto e1 = x1.begin(); e1 != x1.end(); ++e1) {
    auto i = (e1.index * x2_sz);
    nat_t hidx = i / x1.size();
    MatSetValue(M, i, i, 
      τ * hx1[hidx] - β * hx1[hidx] * bsx2[0], ADD_VALUES);
    MatSetValue(M, i + 1, i,      
                  -β * hx1[hidx] * bsx2[1], ADD_VALUES);
    MatSetValue(M, i + 2, i, 
                  -β * hx1[hidx] * bsx2[2], ADD_VALUES);
  }
}

void sbp_sat::x2::write_boundary_x2_right(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bsx2,
  vector_t     const &hx1, 
  real_t       const  β,
  real_t       const  τ,
  int          const  order) {

  auto size = (x1.size() * x2.size()) - 1;
  auto x2_sz = x2.size();

  for (auto e1 = x1.begin(); e1 != x1.end(); ++e1) {
    auto i = (e1.index * x2_sz);

    std::cout << size - i << " " << size - i << std::endl;

    nat_t hidx = i / x1.size();

    MatSetValue(M, size - i, size - i, 
       τ * hx1[hidx] - β * hx1[hidx] * bsx2[0], ADD_VALUES);
    MatSetValue(M, size - i - 1, size - i ,       
                   -β * hx1[hidx] * bsx2[1], ADD_VALUES);
    MatSetValue(M, size - i - 2, size - i, 
                   -β * hx1[hidx] * bsx2[2], ADD_VALUES);
  }
}




