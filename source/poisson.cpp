#include "poisson.h"

void poisson1d::analytical(
  vector_t       &res,
  nat_t  const  size, 
  real_t const  left, 
  real_t const  right) {

  auto h = (right - left) / static_cast<real_t>(size - 1);
  for (std::size_t i = 0; i < size; ++i) 
    res.push_back(std::sin(h * π * i));
}

void poisson1d::analytical_blocked(
  vector_t       &res,
  nat_t    const  size, 
  nat_t    const  local_problems,
  real_t   const  left, 
  real_t   const  right) {

  auto domain_range = range_t(left, right, size);  

  auto local_domain_ranges = std::vector<range_t>();

  auto local_domain_size = (size + local_problems - 1) / local_problems;
  
  /* NOTE: The number of volume points differs from the number of
           matrix points because the interface points are used in 
           both problems. */

  auto volume_points = (local_domain_size - 1) *local_problems+1;

  auto spacing = to_real_t(right - left) / (volume_points - 1);
  for (nat_t i = 0; i != local_problems; ++i) {
    auto begin = volume_points / local_problems * i * spacing;
    auto end   = volume_points / local_problems * (i + 1) * spacing;
    local_domain_ranges.push_back(
      range_t(begin, end, local_domain_size));
  }

  for (std::size_t i = 0; i < local_domain_ranges.size(); ++i) {
    auto local = local_domain_ranges[i];
    for (auto &e : local) {  
      res.push_back(std::sin(π * e));
    }
  }
}

void poisson1d::petsc_problem(
  vector_t         &result,
  nat_t      const  size, 
  nat_t      const  local_problems,
  domain_t   const  domain, 
  boundary_t const  boundary) {
        
  long double β = 1.;
  long double σ1 = -40.;
  long double σ2 = 1;
  long double ϵ = 1.;

  // Unpack tuples 
  auto [left, right] = domain;
  auto [left_data,  left_order, 
        right_data, right_order] = boundary;

  auto domain_range = range_t(left, right, size);  

  auto local_domain_ranges = std::vector<range_t>();

  auto local_domain_size = (size + local_problems - 1) / local_problems;
  
  /* NOTE: The number of volume points differs from the number of
           matrix points because the interface points are used in 
           both problems. */

  auto volume_points = (local_domain_size - 1) *local_problems+1;
  auto matrix_points = local_domain_size * local_problems;

  if (volume_points != size) 
    std::cout << "NOTE: The problem was specified with " 
      << size << " points, but this does not evenly divide " 
      << "the number of local problems (" << local_problems
      <<  "), so it has been adjusted to " << volume_points
      << "." << std::endl;

  auto spacing = to_real_t(right - left) / (volume_points - 1);
  auto spacing_square = spacing * spacing;
  for (nat_t i = 0; i != local_problems; ++i) {
    auto begin = volume_points / local_problems * i * spacing;
    auto end   = volume_points / local_problems * (i + 1) * spacing;
    local_domain_ranges.push_back(
      range_t(begin, end, local_domain_size));
  }

  auto h = numerical::operators::H(local_domain_size, 2, 0, 
    (right - left) / local_problems);
  auto hi = numerical::operators::H_inverse(local_domain_size, 2, 0, 
    (right - left) / local_problems);
  vector_t bs = {(3./2.) / spacing, -2. / spacing, (1./2.) / spacing};
      
  /* Create a petsc object for the A matrix. */
  Mat A;
  MatCreateSeqAIJ(PETSC_COMM_SELF, matrix_points, matrix_points, 
      10, nullptr, &A);

  
  /* Set the second-order SBP D2 in A using D2 sized for 
     the local problem. */
  write_d2(A, local_domain_ranges, local_domain_size, spacing_square);

  write_boundaries(A, local_domain_ranges, hi, bs, matrix_points, 
    local_domain_size, β, σ1, σ2);

  write_fluxs(A, local_domain_ranges, hi, bs, local_domain_size, β, 
    σ1, ϵ);
  
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  Vec b;
  VecCreateSeq(PETSC_COMM_SELF, matrix_points, &b);

  for (std::size_t i = 0; i < local_domain_ranges.size(); ++i) {
    auto local = local_domain_ranges[i];
    auto local_offset = i * local_domain_size;
    for (auto it = local.begin(); it != local.end(); ++it) {  
      auto index = local_offset + it.index;
      VecSetValue(b, index, -π * π * std::sin(*it * π), ADD_VALUES);
    }
  }

  VecSetValue(b, matrix_points - 1, 
    σ2 * right_data * hi[local_domain_size - 1], ADD_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  Vec x;
  VecCreateSeq(PETSC_COMM_SELF, matrix_points, &x);

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPSetUp(ksp);
  KSPSolve(ksp, b, x);
  
  result.resize(matrix_points);

  double *temp;
  VecGetArray(x, &temp);
  for (std::size_t i = 0; i < matrix_points; ++i) 
    result[i] = static_cast<long double>(temp[i]);

  MatDestroy(&A);
  VecDestroy(&b);
  VecDestroy(&x);
}

void poisson1d::petsc_hybridized_problem(
  vector_t         &result,
  nat_t      const  size, 
  nat_t      const  local_problems,
  domain_t   const  domain, 
  boundary_t const  boundary) {

  long double β = 1.;
  long double τ = -40.;
  long double δ_f = 0.0;

  // Unpack tuples 
  auto [left, right] = domain;
  auto [left_data,  left_order, 
        right_data, right_order] = boundary;

  auto interfaces = local_problems - 1;

  auto domain_range = range_t(left, right, size);  

  auto local_domain_ranges = std::vector<range_t>();

  auto local_domain_size = (size + local_problems - 1) / local_problems;

  auto volume_points = (local_domain_size - 1) *local_problems+1;

  if (volume_points != size) 
    std::cout << "NOTE: The problem was specified with " 
      << size << " points, but this does not evenly divide " 
      << "the number of local problems (" << local_problems
      <<  "), so it has been adjusted to " << volume_points
      << "." << std::endl;

  auto spacing = to_real_t(right - left) / (volume_points - 1);
  auto spacing_square = spacing * spacing;
  for (nat_t i = 0; i != local_problems; ++i) {
    auto begin = volume_points / local_problems * i * spacing;
    auto end   = volume_points / local_problems * (i + 1) * spacing;
    local_domain_ranges.push_back(
      range_t(begin, end, local_domain_size));
  }

  auto h = numerical::operators::H(local_domain_size, 2, 0, 
    (right - left) / local_problems);
  auto hi = numerical::operators::H_inverse(local_domain_size, 2, 0, 
    (right - left) / local_problems);
  vector_t bs = {(3./2.) / spacing, -2. / spacing, (1./2.) / spacing};
 
  /*  Create petsc objects for the M block diagonal matrix. 

          [ m_0           ]
      M = [     ...       ]
          [         m_n-1 ]  */

  auto M = std::vector<petsc_matrix>(local_problems); 
  for (auto &m : M) {
    MatCreateSeqAIJ(PETSC_COMM_SELF, local_domain_size, 
      local_domain_size, 4, nullptr, &m);
  }

  /* Initialize each block with a skew diagonal operator. */
  for (auto &m : M) {
    write_d2_h1(m, h, local_domain_ranges, local_domain_size, 
      spacing_square * local_problems);
  }

  /* Write boundary data on each block, first and final blocks potentially
     have unique left and right-hand boundaries accordingly. In this case 
     the final block has a second-order right-hand boundary. */

  write_upper_first_order_boundary(M[0], bs, β, τ);
  write_lower_second_order_boundary(M[local_problems - 1], bs, 
    local_domain_size, β, τ);
  for (nat_t i = 1; i != local_problems; ++i) {
    write_lower_first_order_boundary(M[i - 1], bs, local_domain_size, 
      β, τ);
    write_upper_first_order_boundary(M[i], bs, β, τ);
    
  }
  
  /* Finalize each block (petsc book-keeping). */
  for (auto &m : M) {
    MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
  }

  auto gbar = std::vector<petsc_vector>(local_problems); 
  for (auto &gb : gbar) {
    VecCreateSeq(PETSC_COMM_SELF, local_domain_size, &gb);
  }

  // Set and finalize the g_bar vectors. 
  for (std::size_t i = 0; i < gbar.size(); ++i) {
    auto local = local_domain_ranges[i];
    for (auto it = local.begin(); it != local.end(); ++it) {  
      VecSetValue(gbar[i], it.index, 
        -π * π * std::sin(*it * π) * h[it.index], ADD_VALUES);
    }
    
  }

  VecSetValue(gbar[0], 0, τ * left_data, ADD_VALUES);
  VecSetValue(gbar[0], 0, bs[0] * left_data, ADD_VALUES);
  VecSetValue(gbar[0], 1, bs[1] * left_data, ADD_VALUES);
  VecSetValue(gbar[0], 2, bs[2] * left_data, ADD_VALUES);

  VecSetValue(gbar[local_problems - 1], local_domain_size - 1, 
    right_data, ADD_VALUES);
  VecSetValue(gbar[local_problems - 1], local_domain_size - 1, 
    (1 / τ) * bs[0] * right_data, ADD_VALUES);
  VecSetValue(gbar[local_problems - 1], local_domain_size - 2, 
    (1 / τ) * bs[1] * right_data, ADD_VALUES);
  VecSetValue(gbar[local_problems - 1], local_domain_size - 3, 
    (1 / τ) * bs[2] * right_data, ADD_VALUES);

  for (std::size_t i = 0; i < gbar.size(); ++i) {
    VecAssemblyBegin(gbar[i]);
    VecAssemblyEnd(gbar[i]);
  }

  // F and its transpose are matrices, they can be sliced into one
  // of three types of vectors by the size of the local problem: 
  // left data only, right data only, and empty vectors. Here, F 
  // contains just these simpler vectors.   

  /*  F Transpose:
         right data          left data              empty
      [ 0 ... BS, τ * BS ][ τ * BS, BS ... 0 ][ 0, 0 ... 0, 0, 0 ]
           left data           right data             empty
      [ 0, 0 ... 0, 0, 0 ][ 0 ... BS, τ * BS ][ τ * BS, BS ... 0 ] */

  petsc_vector Fl = {};
  petsc_vector Fr = {};
  VecCreateSeq(PETSC_COMM_SELF, local_domain_size, &Fl);
  VecCreateSeq(PETSC_COMM_SELF, local_domain_size, &Fr);
    
  VecSetValue(Fl, 0, τ + bs[0], ADD_VALUES);
  VecSetValue(Fl, 1,     bs[1], ADD_VALUES);
  VecSetValue(Fl, 2,     bs[2], ADD_VALUES);
  VecSetValue(Fr, local_domain_size - 1, τ + bs[0], ADD_VALUES);
  VecSetValue(Fr, local_domain_size - 2,     bs[1], ADD_VALUES);
  VecSetValue(Fr, local_domain_size - 3,     bs[2], ADD_VALUES);

  VecAssemblyBegin(Fl);
  VecAssemblyEnd(Fl);
  VecAssemblyBegin(Fr);
  VecAssemblyEnd(Fr);

  auto temps = std::vector<petsc_vector>(local_problems * 2); 
  for (auto &t : temps) {
    VecCreateSeq(PETSC_COMM_SELF, local_domain_size, &t);
    VecAssemblyBegin(t);
    VecAssemblyEnd(t);
  }

  auto M_solvers = std::vector<KSP>(local_problems);
  for (std::size_t i = 0; i < local_problems; ++i) { 
    KSPCreate(PETSC_COMM_WORLD, &M_solvers[i]);
    KSPSetOperators(M_solvers[i], M[i], M[i]);
    KSPSetUp(M_solvers[i]);
  }

  for (std::size_t i = 0; i < local_problems; ++i) 
      KSPSolve(M_solvers[i], gbar[i], temps[i]);

  auto gdel = std::vector<double>(interfaces);
  auto temp_scalars = std::vector<double>(local_problems * 2);
  for (std::size_t i = 0; i < local_problems; ++i) {
    VecTDot(Fl, temps[i], &temp_scalars[i * 2]);
    VecTDot(Fr, temps[i], &temp_scalars[i * 2 + 1]);
  }

  /* Sum temp_scalars (ie incomplete gδ values) into gδ. */
  for (std::size_t i = 0; i < interfaces; ++i) {
    gdel[i] = temp_scalars[(i * 2) + 1] + temp_scalars[(i * 2) + 2];
  }

  // TODO: Review this.
  for (auto &e : gdel) e = (interfaces * h[4] * δ_f) - e;

  for (std::size_t i = 0; i < local_problems; ++i) {
    KSPSolve(M_solvers[i], Fl, temps[i * 2]);
    KSPSolve(M_solvers[i], Fr, temps[i * 2 + 1]);
  }

  petsc_matrix λ_denominator = {};

  MatCreateSeqAIJ(PETSC_COMM_SELF, interfaces, interfaces, 2, nullptr, 
    &λ_denominator);
  
  for (auto &e: temp_scalars) e = 0.;


  /*  M \ F => 0 [ 1 Ø ]  ie. { 0, 1, 2, 3, 4, 5 }
                 [ 2 3 ]
                 [ Ø 4 ] 5      
                            =>  [ 1R + 2L + ØØ  ØR + 3L + 4Ø ] [0,0 0,1]
      Ft    =>   [ R L Ø ]      [ 1Ø + 2R + ØL  ØØ + 3R + 4L ] [1,0 1,1]
                 [ Ø R L ]
                                ^^^ This will always be one or two 
                                non-zero values. */ 


  for (std::size_t i = 0; i < interfaces; ++i) {
      
    for (std::size_t j = 0; j < interfaces; ++j) {

      if (i == j) {
        double a, b = 0.;
        VecTDot(Fr, temps[i * 2 + 1], &a);
        VecTDot(Fl, temps[i * 2 + 2], &b);
        MatSetValue(λ_denominator, i, j, -a, ADD_VALUES);
        MatSetValue(λ_denominator, i, j, -b, ADD_VALUES);
      } 
      else if (i - 1 == j) {
        double a = 0.;
        VecTDot(Fr, temps[i * 2], &a);
        MatSetValue(λ_denominator, i, j, -a, ADD_VALUES);
      }
      else if (i + 1 == j) {
        double a = 0.;
        VecTDot(Fl, temps[i * 2 + 3], &a);
        MatSetValue(λ_denominator, i, j, -a, ADD_VALUES);
      }
    }
    MatSetValue(λ_denominator, i, i, 2 * τ, ADD_VALUES);   
  }

  MatAssemblyBegin(λ_denominator, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(λ_denominator, MAT_FINAL_ASSEMBLY);

  petsc_vector λ_numerator = {};
  VecCreateSeq(PETSC_COMM_SELF, interfaces, &λ_numerator);
  for (std::size_t i = 0; i < interfaces; ++i) {
    VecSetValue(λ_numerator, i, gdel[i], ADD_VALUES);
  }

  VecAssemblyBegin(λ_numerator);
  VecAssemblyEnd(λ_numerator);

  // MatView(λ_denominator, PETSC_VIEWER_STDOUT_WORLD);

  petsc_vector λ = {};
  VecCreateSeq(PETSC_COMM_SELF, interfaces, &λ);
  VecAssemblyBegin(λ);
  VecAssemblyEnd(λ);

  KSP λ_solver = {};
  KSPCreate(PETSC_COMM_WORLD, &λ_solver);
  KSPSetOperators(λ_solver, λ_denominator, λ_denominator);
  KSPSetUp(λ_solver);

  // VecView(λ_numerator, PETSC_VIEWER_STDOUT_WORLD);
  
  KSPSolve(λ_solver, λ_numerator, λ);

  double const *λ_data;
  VecGetArrayRead(λ, &λ_data);

  // rhs = (g_bar - F*lambda)
  for (std::size_t i = 0; i < local_problems; ++i) {
    VecCopy(Fl, temps[i * 2]);
    VecCopy(Fr, temps[i * 2 + 1]);
  }

  for (std::size_t i = 0; i < local_problems - 1; ++i) {
    /* NOTE: AXPY => y = a x + y */
    VecAXPY(gbar[i],     -λ_data[i], temps[i * 2 + 1]);
    VecAXPY(gbar[i + 1], -λ_data[i], temps[i * 2 + 2]);
  }

  auto u = std::vector<petsc_vector>(local_problems);
  for (auto &e : u) {
    VecCreateSeq(PETSC_COMM_SELF, local_domain_size, &e);
    VecAssemblyBegin(e);
    VecAssemblyEnd(e);
  }

  for (std::size_t i = 0; i < local_problems; ++i) {
    KSPSolve(M_solvers[i], gbar[i], u[i]);
  }

  // VecView(u[0], PETSC_VIEWER_STDOUT_WORLD);
  // VecView(u[1], PETSC_VIEWER_STDOUT_WORLD);
  // VecView(u[2], PETSC_VIEWER_STDOUT_WORLD);

  for (auto &e : u) {
    double const *u_data;
    VecGetArrayRead(e, &u_data);
    for (std::size_t i = 0; i < local_domain_size; ++i) 
      result.push_back(u_data[i]);
  }

  for (auto &e : M)         MatDestroy(&e);
  for (auto &e : gbar)      VecDestroy(&e);
  for (auto &e : M_solvers) KSPDestroy(&e);
  VecDestroy(&Fl);
  VecDestroy(&Fr);
  for (auto &e : temps)     VecDestroy(&e);
  VecDestroy(&λ_numerator);
  MatDestroy(&λ_denominator);
  KSPDestroy(&λ_solver);
  VecDestroy(&λ);
  for (auto &e : u)         VecDestroy(&e);
}

void poisson1d::write_d2(
  petsc_matrix               &A, 
  std::vector<range_t> const &local_domain_ranges,
  long double          const  local_domain_size,
  long double          const  spacing_square) {
  
  for (std::size_t i = 0; i < local_domain_ranges.size(); ++i) {
    auto local = local_domain_ranges[i];
    auto local_offset = i * local_domain_size;
    
    /* Initialize the first local skew row. */
    MatSetValue(A, local_offset, local_offset,  
       1. / spacing_square, ADD_VALUES);
    MatSetValue(A, local_offset, local_offset + 1, 
      -2. / spacing_square, ADD_VALUES);
    MatSetValue(A, local_offset, local_offset + 2, 
       1. / spacing_square, ADD_VALUES); 

    /* Initialize the final local skew row. */
    auto local_last_row = (i + 1) * local_domain_size;
    MatSetValue(A, local_last_row - 1, local_last_row - 3,  
       1. / spacing_square, ADD_VALUES);
    MatSetValue(A, local_last_row - 1, local_last_row - 2, 
      -2. / spacing_square, ADD_VALUES);
    MatSetValue(A, local_last_row - 1, local_last_row - 1, 
       1. / spacing_square, ADD_VALUES); 

    /* Initialize the interior local diagonal rows. */
    for (auto it = local.begin() + 1; it != local.end() - 1; ++it) {  
      auto global_row_index    = local_offset + it.index;
      auto global_column_index = local_offset + it.index - 1;

      MatSetValue(A, global_row_index, global_column_index,  
         1. / spacing_square, ADD_VALUES);
      MatSetValue(A, global_row_index, global_column_index + 1, 
        -2. / spacing_square, ADD_VALUES);
      MatSetValue(A, global_row_index, global_column_index + 2, 
         1. / spacing_square, ADD_VALUES);  
    }
  }
}

void poisson1d::write_d2_h1(
  petsc_matrix                   &M, 
  std::vector<long double> const &h, 
  std::vector<range_t>     const &local_domain_ranges,
  long double              const  local_domain_size,
  long double              const  spacing_square) {
  
  for (std::size_t i = 0; i < local_domain_ranges.size(); ++i) {
    auto local = local_domain_ranges[i];
    
    /* Initialize the first local skew row. */
    MatSetValue(M, 0, 0,  1. / spacing_square * h[0], ADD_VALUES);
    MatSetValue(M, 0, 1, -2. / spacing_square * h[0], ADD_VALUES);
    MatSetValue(M, 0, 2,  1. / spacing_square * h[0], ADD_VALUES); 

    /* Initialize the final local skew row. */
    auto n = local_domain_size - 1;
    MatSetValue(M, n, n - 2,  1. / spacing_square * h[n], ADD_VALUES);
    MatSetValue(M, n, n - 1, -2. / spacing_square * h[n], ADD_VALUES);
    MatSetValue(M, n, n,      1. / spacing_square * h[n], ADD_VALUES); 

    /* Initialize the interior local diagonal rows. */
    for (auto it = local.begin() + 1; it != local.end() - 1; ++it) {  
      auto ij = it.index;
      MatSetValue(M, ij, ij - 1,  1./spacing_square*h[ij], ADD_VALUES);
      MatSetValue(M, ij, ij,     -2./spacing_square*h[ij], ADD_VALUES);
      MatSetValue(M, ij, ij + 1,  1./spacing_square*h[ij], ADD_VALUES);  
    }
  }
}

void poisson1d::write_boundaries(
    petsc_matrix               &A, 
    std::vector<range_t> const &local_domain_ranges,
    vector_t             const &hi, 
    vector_t             const &bs, 
    long double          const  matrix_points,
    long double          const  local_domain_size,
    long double          const  β,
    long double          const  σ1,
    long double          const  σ2) {

  /* Set the second-order SAT terms for the left hand first-order
       boundary: β * HI1 * BS' * e0 * e0' */
  MatSetValue(A, 0, 0, β * hi[0] * bs[0], ADD_VALUES);
  MatSetValue(A, 1, 0, β * hi[1] * bs[1], ADD_VALUES);
  MatSetValue(A, 2, 0, β * hi[2] * bs[2], ADD_VALUES);

  /* Set the boundary data for the left hand first-order boundary. */
  MatSetValue(A, 0, 0, σ1 * hi[0], ADD_VALUES);

  MatSetValue(A, matrix_points - 1, matrix_points - 3, 
    σ2 * hi[local_domain_size - 1] * bs[2], ADD_VALUES);
  MatSetValue(A, matrix_points - 1, matrix_points - 2, 
    σ2 * hi[local_domain_size - 1] * bs[1], ADD_VALUES);
  MatSetValue(A, matrix_points - 1, matrix_points - 1, 
    σ2 * hi[local_domain_size - 1] * bs[0], ADD_VALUES);

}

void poisson1d::write_fluxs(
    petsc_matrix               &A, 
    std::vector<range_t> const &local_domain_ranges,
    vector_t             const &hi, 
    vector_t             const &bs, 
    long double          const  local_domain_size,
    long double          const  β,
    long double          const  σ1,
    long double          const  ϵ) {

  for (std::size_t i = 0; i < local_domain_ranges.size() - 1; ++i) {
    auto local = local_domain_ranges[i];
    auto local_offset = (i + 1) * local_domain_size;
    auto n = local_domain_size;

    // β * HI1 * BS' * en * en' 
    MatSetValue(A, local_offset - 3, local_offset - 1, 
      β * hi[n - 3] * bs[2], ADD_VALUES);
    MatSetValue(A, local_offset - 2, local_offset - 1, 
      β * hi[n - 2] * bs[1], ADD_VALUES);
    MatSetValue(A, local_offset - 1, local_offset - 1, 
      β * hi[n - 1] * bs[0], ADD_VALUES);

    // ϵ * HI1 * BS  * en * en' 
    MatSetValue(A, local_offset - 1, local_offset - 3, 
      ϵ * hi[n - 1] * bs[2], ADD_VALUES);
    MatSetValue(A, local_offset - 1, local_offset - 2, 
      ϵ * hi[n - 1] * bs[1], ADD_VALUES);
    MatSetValue(A, local_offset - 1, local_offset - 1, 
      ϵ * hi[n - 1] * bs[0], ADD_VALUES);

    // σ₁ * HI1 * en * en' 
    MatSetValue(A, local_offset - 1, local_offset - 1, 
      σ1 * hi[n - 1], ADD_VALUES);

    // β * HI1 * BS' *e0 * e0' 
    MatSetValue(A, local_offset, local_offset, 
      β * hi[0] * bs[0], ADD_VALUES);
    MatSetValue(A, local_offset + 1, local_offset, 
      β * hi[1] * bs[1], ADD_VALUES);
    MatSetValue(A, local_offset + 2, local_offset, 
      β * hi[2] * bs[2], ADD_VALUES);

    // ϵ * HI1 * BS * e0 * e0 
    MatSetValue(A, local_offset, local_offset, 
      ϵ * hi[0] * bs[0], ADD_VALUES);
    MatSetValue(A, local_offset, local_offset + 1, 
      ϵ * hi[0] * bs[1], ADD_VALUES);
    MatSetValue(A, local_offset, local_offset + 2, 
      ϵ * hi[0] * bs[2], ADD_VALUES);

    // σ₁ * HI1 * e0 * e0' 
    MatSetValue(A, local_offset, local_offset, 
      σ1 * hi[0], ADD_VALUES);

    // -σ₁ * HI1 * en * e0' 
    MatSetValue(A, local_offset - 1, local_offset, 
      -σ1 * hi[n - 1], ADD_VALUES);

    // ϵ * HI1 * en * e0' * BS 
     MatSetValue(A, local_offset - 1, local_offset, 
      ϵ * hi[n - 1] * bs[0], ADD_VALUES);
    MatSetValue(A, local_offset - 1, local_offset + 1, 
      ϵ * hi[n - 1] * bs[1], ADD_VALUES);
    MatSetValue(A, local_offset - 1, local_offset + 2, 
      ϵ * hi[n - 1] * bs[2], ADD_VALUES);

    // -β * HI1 * BS' * en * e0' 
    MatSetValue(A, local_offset - 3, local_offset, 
      -β * hi[n - 3] * bs[2], ADD_VALUES);
    MatSetValue(A, local_offset - 2, local_offset, 
      -β * hi[n - 2] * bs[1], ADD_VALUES);
    MatSetValue(A, local_offset - 1, local_offset, 
      -β * hi[n - 1] * bs[0], ADD_VALUES);

    // -σ₁ * HI1 * e0 * en' 
    MatSetValue(A, local_offset, local_offset - 1, 
      -σ1 * hi[0], ADD_VALUES);

    // ϵ * HI1 * e0 * en' * BS 
     MatSetValue(A, local_offset, local_offset - 1, 
      ϵ * hi[0] * bs[0], ADD_VALUES);
    MatSetValue(A, local_offset, local_offset - 2, 
      ϵ * hi[0] * bs[1], ADD_VALUES);
    MatSetValue(A, local_offset, local_offset - 3, 
      ϵ * hi[0] * bs[2], ADD_VALUES);

    // -β * HI1 * BS' * e0 * en' 
    MatSetValue(A, local_offset + 2, local_offset - 1, 
      -β * hi[2] * bs[2], ADD_VALUES);
    MatSetValue(A, local_offset + 1, local_offset - 1, 
      -β * hi[1] * bs[1], ADD_VALUES);
    MatSetValue(A, local_offset, local_offset - 1,
      -β * hi[0] * bs[0], ADD_VALUES);      
  }
}


void poisson1d::write_upper_first_order_boundary(
  petsc_matrix       &M, 
  vector_t     const &bs, 
  long double  const  β,
  long double  const  τ) {
      
  MatSetValue(M, 0, 0, τ, ADD_VALUES);
  MatSetValue(M, 0, 0, β * bs[0], ADD_VALUES);
  MatSetValue(M, 1, 0, β * bs[1], ADD_VALUES);
  MatSetValue(M, 2, 0, β * bs[2], ADD_VALUES);
}

void poisson1d::write_lower_first_order_boundary(
  petsc_matrix       &M, 
  vector_t     const &bs,
  long double  const  n, 
  long double  const  β,
  long double  const  τ) {
    
  MatSetValue(M, n - 1, n - 1, τ, ADD_VALUES);
  MatSetValue(M, n - 1, n - 1, β * bs[0], ADD_VALUES);
  MatSetValue(M, n - 2, n - 1, β * bs[1], ADD_VALUES);
  MatSetValue(M, n - 3, n - 1, β * bs[2], ADD_VALUES);
}

void poisson1d::write_lower_second_order_boundary(
  petsc_matrix       &M, 
  vector_t     const &bs,
  long double  const  n, 
  long double  const  β,
  long double  const  τ) {
    
  MatSetValue(M, n - 1, n - 1, β * bs[0], ADD_VALUES);
  MatSetValue(M, n - 1, n - 2, β * bs[1], ADD_VALUES);
  MatSetValue(M, n - 1, n - 3, β * bs[2], ADD_VALUES);

  MatSetValue(M, n - 1, n - 1, (1 / τ) * bs[0] * bs[0], ADD_VALUES);
  MatSetValue(M, n - 1, n - 2, (1 / τ) * bs[1] * bs[0], ADD_VALUES);
  MatSetValue(M, n - 1, n - 3, (1 / τ) * bs[2] * bs[0], ADD_VALUES);

  MatSetValue(M, n - 2, n - 1, (1 / τ) * bs[0] * bs[1], ADD_VALUES);
  MatSetValue(M, n - 2, n - 2, (1 / τ) * bs[1] * bs[1], ADD_VALUES);
  MatSetValue(M, n - 2, n - 3, (1 / τ) * bs[2] * bs[1], ADD_VALUES);

  MatSetValue(M, n - 3, n - 1, (1 / τ) * bs[0] * bs[2], ADD_VALUES);
  MatSetValue(M, n - 3, n - 2, (1 / τ) * bs[1] * bs[2], ADD_VALUES);
  MatSetValue(M, n - 3, n - 3, (1 / τ) * bs[2] * bs[2], ADD_VALUES);
}
