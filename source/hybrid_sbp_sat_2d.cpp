#include "hybrid_sbp_sat_2d.h"

/* [ M  F ][ u ] = [ g   ]
     [ FT D ][ λ ]   [ g_λ ]
  
    [ B   T ][ u ] = [ bu ]
    [ TT  D ][ λ ]   [ bλ ]
  */

/* Baseline Poisson solution. Solves A u = b given a 2D volume point 
   dimension and boundary conditions. */
void sbp_sat::x2::
petsc_poisson(sbp_sat::real_v             &result,
              sbp_sat::domain_v     const &domain, 
              sbp_sat::block_t      const &block, 
              sbp_sat::boundary_tx2 const &boundaries) {
 

  using namespace sbp_sat;

  auto A = matrix<fw> {};

  auto w = neumann;
  auto e = neumann;
  auto n = dirichlet;
  auto s = dirichlet;

  x2::write_m(A, block, block, {w, e, n, s});

  finalize<fw>(A);

  destroy<fw>(A);

}

void sbp_sat::x2::
petsc_hybridized_poisson(sbp_sat::real_v             &result,
                         sbp_sat::domain_v     const &domain, 
                         sbp_sat::block_v      const &blocks, 
                         sbp_sat::boundary_vx2 const &boundaries) {



  /* [ M  F ][ u ] = [ g   ]
     [ FT D ][ λ ]   [ g_λ ]
  
    [ B   T ][ u ] = [ bu ]
    [ TT  D ][ λ ]   [ bλ ]
  */

  logging::out << "begin." << std::endl;
   
  using namespace sbp_sat;

  nat_t local_problems = blocks.size() * blocks.size();

  // only for square hybrid systems 
  const std::size_t n_blocks = blocks.size() * blocks.size();

  /*  Create petsc objects for the block diagonal matrix, M. 

          [ m_0           ] 
      M = [     ...       ] in ( M  F )  
          [         m_n-1 ]    ( Ft D )   

      From the segmented volume point blocks

      | b0  | ... | ...  |
      |-----|-----|------|
      | ... | ... | ...  |
      |-----|-----|------|
      | ... | ... | bn-1 |  */

  auto gw = [](real_t x, real_t y){return std::sin(π * x + π * y);};
  auto ge = [](real_t x, real_t y){return std::sin(π * x + π * y);};
  auto gs = [](real_t x, real_t y){return -π * std::cos(π * x + π * y);};
  auto gn = [](real_t x, real_t y){return π * std::cos(π * x + π * y);};


  // Add a 2 here because u_xx == u_yy
  auto source = [](real_t x, real_t y) {
    return 2. * (-π * π * sin(π * x + π * y));};


  std::vector<petsc_matrix> sources;

  std::vector<range_t> grids = {
    range_t(0, 1./3., 4), 
    range_t(1./3., 2./3., 4), 
    range_t(2./3., 1., 4)};

  compute_sources(sources, grids, source);

  std::vector<std::vector<std::size_t>> interfaces = {
    {0, 1, 0, 3, 0, 0, 0, 0, 0}, 
    {0, 0, 2, 0, 4, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 5, 0, 0, 0}, 
    {0, 0, 0, 0, 6, 0, 8, 0, 0}, 
    {0, 0, 0, 0, 0, 7, 0, 9, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 10}, 
    {0, 0, 0, 0, 0, 0, 0, 11, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 12}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0}};
/*
   {{0, 1, 0, 7, 0, 0, 0, 0, 0}, 
    {0, 0, 2, 0, 8, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 9, 0, 0, 0}, 
    {0, 0, 0, 0, 3, 0, 10, 0, 0}, 
    {0, 0, 0, 0, 0, 4, 0, 11, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 12}, 
    {0, 0, 0, 0, 0, 0, 0, 5, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 6}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0}};



    {0, 1, 0, 3, 0, 0, 0, 0, 0}, 
    {0, 0, 2, 0, 4, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 5, 0, 0, 0}, 
    {0, 0, 0, 0, 6, 0, 7, 0, 0}, 
    {0, 0, 0, 0, 0, 8, 0, 9, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 10}, 
    {0, 0, 0, 0, 0, 0, 0, 11, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 12}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0}};

  1   0   4.0   0    0    0    0    0    0    0    ⋅    ⋅ 
   2  1.0   ⋅   4.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅   2.0   ⋅    ⋅   4.0  0.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅    ⋅   3.0   ⋅   0.0  1.0   ⋅   4.0   ⋅    ⋅    ⋅    ⋅ 
  ⋅    ⋅    ⋅   3.0  0.0  2.0  1.0   ⋅   4.0   ⋅    ⋅    ⋅ 
  ⋅    ⋅    ⋅    ⋅   3.0  0.0  2.0   ⋅    ⋅   4.0  0.0   ⋅ 
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   3.0   ⋅   0.0  1.0   ⋅ 
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   3.0  0.0  2.0  1.0
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   3.0  0.0  2.0


*/

  std::size_t n_interfaces = 0;
  for (auto &row : interfaces) {
    for (auto &element : row) {
      if (element != 0) {
        n_interfaces += 1;
      }
    }
  }

  auto [rows, cols] = blocks[0];
  auto spacing = rows.to - rows.from;
  auto span = spacing;
  auto n_points_x = rows.size();
  // components class contains most of the sbp-sat component matrices
  // needed to set up an sbp-sat problem.
  auto sbp = components{4, span};

  sbp.τ = 42.; // hard code these coeffs for now. 
  sbp.β = 1.;

  sbp.n_blocks = n_blocks;         // additional non sbp-sat info. but  
  sbp.n_interfaces = n_interfaces; // useful to have along. 

  

  // Compute flux components used by b. 
  std::vector<petsc_matrix> B;
  compute_B(B, sbp);  

  // Generate the solution at the boundary. This implementation generates 
  // vectors from range functions. 
  vv<petsc_vector> boundary_solution;
  compute_boundary_solution(  
    boundary_solution, grids, {gw, ge, gs, gn}, {0., 1., 0., 1.});

  // Maps blocks with boundary solutions 
  vv<std::size_t> boundary_data_map = {
    {1, 0, 1, 0}, 
    {2, 0, 0, 0}, 
    {3, 0, 0, 1},
    {0, 0, 2, 0}, 
    {0, 0, 0, 0}, 
    {0, 0, 0, 2}, 
    {0, 1, 3, 0}, 
    {0, 2, 0, 0}, 
    {0, 3, 0, 3}};

  vv<std::size_t> boundary_order_map = {
    {1, 0, 2, 0}, 
    {1, 0, 0, 0}, 
    {1, 0, 0, 2},
    {0, 0, 2, 0}, 
    {0, 0, 0, 0}, 
    {0, 0, 0, 2}, 
    {0, 1, 2, 0}, 
    {0, 1, 0, 0}, 
    {0, 1, 0, 2}};

  // Compute the hybrid system g terms. 
  std::vector<petsc_vector> g;
  compute_g(g, B, boundary_solution, sources, boundary_order_map, 
    boundary_data_map, sbp); 

  // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  // MatView(B[7], PETSC_VIEWER_STDOUT_SELF);

  // 2: λ_rhs    = gd - FT * M \ g 






  vv<std::size_t> F_symbols(n_blocks, // rows 
    std::vector<std::size_t>(n_interfaces, 0)); // columns
  vv<std::size_t> FT_symbols(n_interfaces, // rows
    std::vector<std::size_t>(n_blocks, 0)); // columns

  sbp_sat::x2::initialize_symbols(F_symbols, FT_symbols, interfaces, sbp);

  

  std::cout << "Square hybrid specs:" << std::endl 
    << " | local problem size: " << n_points_x << " x " << n_points_x 
    << std::endl << " | span: " << span 
    << std::endl << " | total blocks: " << n_blocks 
    << std::endl << " | total interfaces: " << n_interfaces 
    << std::endl;


  std::cout << n_points_x << ", " << span << std::endl;

  // n_points = number of columns in 
  //            matrix component. 

  //petsc_vector λb
  //initialize_λb(λb, sbp);


  auto M = std::vector<petsc_matrix>(sbp.n_blocks);

  // hard code for now
  sbp_sat::x2::make_M(M[0], sbp, {1, 1, 2, 1});
  sbp_sat::x2::make_M(M[1], sbp, {1, 1, 1, 1});
  sbp_sat::x2::make_M(M[2], sbp, {1, 1, 1, 2});
  sbp_sat::x2::make_M(M[3], sbp, {1, 1, 2, 1});
  sbp_sat::x2::make_M(M[4], sbp, {1, 1, 1, 1});
  sbp_sat::x2::make_M(M[5], sbp, {1, 1, 1, 2});
  sbp_sat::x2::make_M(M[6], sbp, {1, 1, 2, 1});
  sbp_sat::x2::make_M(M[7], sbp, {1, 1, 1, 1});
  sbp_sat::x2::make_M(M[8], sbp, {1, 1, 1, 2});

  // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  // MatView(M[0], PETSC_VIEWER_STDOUT_SELF);

  // Make direct solvers
  auto solvers = std::vector<KSP>(M.size());
  for (std::size_t i = 0; i != M.size(); ++i) {
    solvers[i] = KSP();
    KSPCreate(PETSC_COMM_SELF, &solvers[i]); 
    KSPSetOperators(solvers[i], M[i], M[i]);
  }

  // Make F components
  // (stored as vectors because we compute M^-1 F by solving every row 
  //  for each component.) 
  auto f = std::vector<std::vector<petsc_vector>>(
    4, std::vector<petsc_vector>(n_blocks));  
  auto F = std::vector<petsc_matrix>(4);  
  make_F(sbp, F, f);

  auto MF = vv<petsc_vector>(4 * sbp.n_blocks, 
            std::vector<petsc_vector>(sbp.n)); 
  initialize_MF(MF, sbp);

  compute_MF(MF, solvers, f);

  // MF {M0 FW, M0 FE, M0 FS, M0 MN, M1 FW, ..., MN FN}

  petsc_matrix D;
  compute_D(D, sbp, interfaces);
  
  
  // Check that decomposed MF: 
  //   - is numerically correct and 
  //   - aligns with the super indices.  
  // print_MF(MF, F_symbols, sbp);
  /*
  
  
  auto FTMF = vv<petsc_matrix>(n_interfaces,
    std::vector<petsc_matrix>(n_interfaces));
  for (std::size_t i = 0; i != n_interfaces; ++i) {
    for (std::size_t j = 0; j != n_interfaces; ++j) {
        MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n, sbp.n, 
          sbp.n, nullptr, &FTMF[i][j]);
    }
  }
  compute_ftmf(FTMF, f, MF, F_symbols, FT_symbols, sbp);
  */

  petsc_matrix λA; 
  initialize_λA(λA, sbp);
  compute_λA(λA, D, f, MF, F_symbols, FT_symbols, sbp);

  logging::out << "compute Mg..." << std::endl;

  std::vector<petsc_vector> Mg;
  initialize_Mg(Mg, sbp);
  compute_Mg(Mg, solvers, g);

  logging::out << "compute ..." << std::endl;

  petsc_vector λb;
  initialize_λb(λb, sbp); 
  // λb = -FT * M \ g 
  compute_λb(λb, F, Mg, FT_symbols, sbp);

  logging::out << "exit..." << std::endl;
  exit(-1);

  destroy_MF(MF);
  destroy_Mg(Mg);
  destroy<fw>(λA);

  for (auto &e : M) destroy<fw>(e);


  // print_FTMF(f, MF, F_symbols, FT_symbols, sbp);


  //MatView(FTMF[0][0], PETSC_VIEWER_STDOUT_SELF);

  // lambda_matrix = compute_dftmf()  

  // 1: λ_matrix =  D - FT * M \ F 
  // 2: λ_rhs    = gd - FT * M \ g 
  // 3: λ = λ_matrix \ λ_rhs

  // compute_mg() 

  logging::out << "end." << std::endl;
}

void sbp_sat::x2::
write_m(
  petsc_matrix       &m,
  block_t      const &block_x1, 
  block_t      const &block_x2, 
  std::array<std::size_t, 4> const bc) {

  auto [rows, c_] = block_x1;
  auto [r_, cols] = block_x2;
  (void) c_; (void) r_; // unused;

  auto x1_spacing = rows.to - rows.from;
  auto x2_spacing = cols.to - cols.from;

  auto x1_sz = rows.size() - 1;
  auto x2_sz = cols.size() - 1;

  auto x1_spacing_square = x1_spacing * x1_spacing;
  auto x2_spacing_square = x2_spacing * x2_spacing;


  auto hx1 = numerical::operators::H(rows.size(), 2, 0, x1_spacing);
  auto hx2 = numerical::operators::H(cols.size(), 2, 0, x2_spacing);

  // hard code for now 
  real_t τ = 42.;
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
  MatCreateSeqAIJ(PETSC_COMM_SELF, rows.size(), rows.size(), 4, 
    nullptr, &d2x1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, cols.size(), cols.size(), 4, 
    nullptr, &d2x2);
  MatCreateSeqAIJ(PETSC_COMM_SELF, rows.size(), rows.size(), 1, 
    nullptr, &h1x1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, cols.size(), cols.size(), 1, 
    nullptr, &h1x2);

  write_d2_h1(d2x1, hx1, rows, x1_spacing_square / (x1_sz * x1_sz), -1.);
  write_d2_h1(d2x2, hx2, cols, x2_spacing_square / (x2_sz * x2_sz), -1.);
  write_h1(h1x1, hx1);
  write_h1(h1x2, hx2);
  finalize<fw>(d2x1);
  finalize<fw>(d2x2);
  finalize<fw>(h1x1);
  finalize<fw>(h1x2);

  MatSeqAIJKron(h1x2, d2x1, MAT_INITIAL_MATRIX, &Ax1);
  MatSeqAIJKron(d2x2, h1x1, MAT_INITIAL_MATRIX, &Ax2);

  finalize<fw>(Ax2);
  finalize<fw>(Ax1);
  
  auto sz = rows.size() * cols.size();
  MatCreate(PETSC_COMM_SELF, &m);
  MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, sz, sz);
  MatSetType(m, MATCOMPOSITE);

  MatCompositeAddMat(m, Ax1);
  MatCompositeAddMat(m, Ax2);

  // Append boundary condition coefficients as new matrices to the 
  // composite matrix, m. 
  add_boundary<x, left> (m, rows, cols, bsx1, hx2, β, τ, bc[0]);
  add_boundary<x, right>(m, rows, cols, bsx1, hx2, β, τ, bc[1]);
  add_boundary<y, left> (m, rows, cols, bsx1, hx2, β, τ, bc[2]);
  add_boundary<y, right>(m, rows, cols, bsx1, hx2, β, τ, bc[3]);
  
  MatCompositeSetMatStructure(m, DIFFERENT_NONZERO_PATTERN);
  MatCompositeMerge(m);
  finalize<fw>(m);

  destroy<fw>(d2x1);
  destroy<fw>(d2x2);
  destroy<fw>(h1x1);
  destroy<fw>(h1x2);
  destroy<fw>(Ax1);
  destroy<fw>(Ax2);
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

void sbp_sat::x2::write_Ls(
  std::vector<sbp_sat::petsc_matrix> &L,
  sbp_sat::boundary_vx2 const &boundaries){
  for (std::size_t i = 0; i != boundaries.size(); ++i) {
    for (std::size_t j = 0; j != boundaries[i].size(); ++j) {
      auto [g_t, _] = boundaries[i][j];
      auto ell = &L[i + j * boundaries.size()];
      std::cout << i + j * boundaries.size() << std::endl;
      MatCreateSeqAIJ(PETSC_COMM_SELF, g_t.size(), 
        g_t.size() * g_t.size(), 4, nullptr, ell);

      if (j == 0) {
        for (std::size_t k = 0; k != g_t.size(); ++k)
          MatSetValue(*ell, k, k,  1., ADD_VALUES);
      }
      else if (j == 1) {
        for (std::size_t k = 0; k != g_t.size(); ++k)
          MatSetValue(*ell, k, g_t.size() * (g_t.size() - 1) + k,  1., ADD_VALUES);
      }
      else if (j == 2) {
        for (std::size_t k = 0; k != g_t.size(); ++k)
          MatSetValue(*ell, k, k * g_t.size(),  1., ADD_VALUES);
      }
      else if (j == 3) {
        for (std::size_t k = 0; k != g_t.size(); ++k)
          MatSetValue(*ell, k, (k + 1) * g_t.size() - 1,  1., ADD_VALUES);
      }
      linalg::finalize<sbp_sat::x2::fw>(*ell);
    }   
  }
}

void write_Lts(
  std::vector<sbp_sat::petsc_matrix> &L,
  std::vector<sbp_sat::petsc_matrix> &Lt) { 

    for (auto i = std::size_t(0); i != L.size(); ++i) {
      MatTranspose(L[i], MAT_INITIAL_MATRIX, &Lt[Lt.size() - i]);
      linalg::finalize<sbp_sat::x2::fw>(Lt[Lt.size() - i]);
    }
}

void sbp_sat::x2::solve(
  KSP &A, std::vector<petsc_vector> &b, std::vector<petsc_vector> &x) {
    for (std::size_t i = 0; i < b.size(); ++i) {
        KSPSolve(A, b[i], x[i]);
    }
}

void sbp_sat::x2::fcompop(
  petsc_matrix &f, 
  petsc_matrix const &l, 
  petsc_matrix const &b,
  petsc_matrix const &h,
  real_t       const τ, 
  real_t       const β) {

  petsc_matrix t;
  MatMatMult(l, b, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &t);
  MatScale(t, β);
  MatAXPY(t, -τ, l, UNKNOWN_NONZERO_PATTERN);
  MatMatMult(t, h, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &f);
  finalize<fw>(f);
  destroy<fw>(t);
}

void sbp_sat::x2::make_F(
  components const &sbp, 
  std::vector<petsc_matrix> &F,
  std::vector<std::vector<petsc_vector>> &f) {

  // (-τ * LN + β * LN* BS_y) * H_x 
  fcompop(F[3], sbp.ln, sbp.bsy, sbp.hx, sbp.τ, sbp.β);
  // (-τ * LS + β * LS* BS_x) * H_x 
  fcompop(F[2], sbp.ls, sbp.bsy, sbp.hx, sbp.τ, sbp.β);
  // (-τ * LE + β * LE* BS_y) * H_y 
  fcompop(F[1], sbp.le, sbp.bsx, sbp.hy, sbp.τ, sbp.β);
  // (-τ * LW + β * LW* BS_y) * H_y 
  fcompop(F[0], sbp.lw, sbp.bsx, sbp.hy, sbp.τ, sbp.β);

  int ncols;
  int const *cols;
  const double *vals;

  f[0] = std::vector<petsc_vector>(sbp.n);
  f[1] = std::vector<petsc_vector>(sbp.n);
  f[2] = std::vector<petsc_vector>(sbp.n);
  f[3] = std::vector<petsc_vector>(sbp.n);

  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);

  for (std::size_t i = 0; i != sbp.n; ++i) {
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &f[0][i]);
    MatGetRow(F[0], i, &ncols, &cols, &vals);
    VecSetValues(f[0][i], ncols, cols, vals, ADD_VALUES);
    finalize<fw>(f[0][i]);

    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &f[1][i]);
    MatGetRow(F[1], i, &ncols, &cols, &vals);
    VecSetValues(f[1][i], ncols, cols, vals, ADD_VALUES);
    finalize<fw>(f[1][i]);
  
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &f[2][i]);
    MatGetRow(F[2], i, &ncols, &cols, &vals);
    VecSetValues(f[2][i], ncols, cols, vals, ADD_VALUES);
    finalize<fw>(f[2][i]);
  
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &f[3][i]);
    MatGetRow(F[3], i, &ncols, &cols, &vals);
    VecSetValues(f[3][i], ncols, cols, vals, ADD_VALUES);
    finalize<fw>(f[3][i]);
  }
}

void sbp_sat::x2::make_M(
  petsc_matrix &M,
  components const &sbp, 
  std::array<std::size_t, 4> const &boundary) {

  petsc_matrix temp1, temp2; 
  MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, nullptr, 
    &temp1);
  MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, nullptr, 
    &temp2);
  finalize<fw>(temp1);
  finalize<fw>(temp2);

  MatAXPY(temp1, 1., sbp.d2x, UNKNOWN_NONZERO_PATTERN);
  MatAXPY(temp1, 1., sbp.d2y, UNKNOWN_NONZERO_PATTERN);
  MatAXPY(temp2, -1., sbp.hl, UNKNOWN_NONZERO_PATTERN);

  MatMatMult(temp2, temp1, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &M);

  make_M_boundary(M, sbp, 1, boundary[0]);
  make_M_boundary(M, sbp, 2, boundary[1]);
  make_M_boundary(M, sbp, 3, boundary[2]);
  make_M_boundary(M, sbp, 4, boundary[3]);

  finalize<fw>(M);
}

void sbp_sat::x2::make_M_boundary(
  petsc_matrix &M,
  components const &sbp, 
  std::size_t const direction,
  std::size_t const boundary) {

  auto L = direction == 1 ? sbp.lw 
         : direction == 2 ? sbp.le
         : direction == 3 ? sbp.ls 
         : sbp.ln; 

  auto H = direction < 3 ? sbp.hy : sbp.hx;
  auto BS = direction < 3 ? sbp.bsx : sbp.bsy;

  if (boundary == dirichlet) {

    petsc_matrix LT, BST, temp1, temp2, temp3, temp4, temp5, temp6, 
      temp7;

    MatTranspose(L, MAT_INITIAL_MATRIX, &LT);
    MatTranspose(BS, MAT_INITIAL_MATRIX, &BST);
    finalize<fw>(LT);
    finalize<fw>(BST);

    // τ*H_y*LW'*LW 
    MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, 
      nullptr, &temp1);
    finalize<fw>(temp1);    
    MatAXPY(temp1, sbp.τ, H, UNKNOWN_NONZERO_PATTERN);    
    MatMatMult(LT, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2);
    finalize<fw>(temp2);    
    MatMatMult(temp1, temp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp3);
    finalize<fw>(temp3);
    MatAXPY(M, 1., temp3, UNKNOWN_NONZERO_PATTERN);

    // -β*H_y*BS_x'*LW'*LW
    MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, 
      nullptr, &temp4);
    finalize<fw>(temp4);
    MatAXPY(temp4, sbp.β, H, UNKNOWN_NONZERO_PATTERN);
    MatMatMult(temp4, BST, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp5);
    finalize<fw>(temp5);
    MatMatMult(temp5, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp6);
    finalize<fw>(temp6);
    MatMatMult(temp6, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp7);
    finalize<fw>(temp7);
    MatAXPY(M, -1., temp7, UNKNOWN_NONZERO_PATTERN);

    // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(temp7, PETSC_VIEWER_STDOUT_SELF);

    destroy<fw>(LT);
    destroy<fw>(BST);
    destroy<fw>(temp1);
    destroy<fw>(temp2);
    destroy<fw>(temp3);
    destroy<fw>(temp4);
    destroy<fw>(temp5);
    destroy<fw>(temp6);
    destroy<fw>(temp7);

  }
  else {

    petsc_matrix LT, BST, temp1, temp2, temp3, temp4, temp5, temp6, 
      temp7, temp8;

    MatTranspose(L, MAT_INITIAL_MATRIX, &LT);
    MatTranspose(BS, MAT_INITIAL_MATRIX, &BST);
    finalize<fw>(LT);
    finalize<fw>(BST);

    // H_x * LN' * LN * BS_y
    MatMatMult(H, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp1);
    finalize<fw>(temp1);
    MatMatMult(temp1, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2);
    finalize<fw>(temp2);
    MatMatMult(temp2, BS, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp3);
    finalize<fw>(temp3);
    MatAXPY(M, 1., temp3, UNKNOWN_NONZERO_PATTERN);

    // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(temp3, PETSC_VIEWER_STDOUT_SELF);

    // -1/τ * H_x * BS_y' * LN' * LN * BS_y 
    MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n * sbp.n, sbp.n * sbp.n, sbp.n, 
      nullptr, &temp4);
    finalize<fw>(temp4);
    MatAXPY(temp4, 1. / sbp.τ, H, UNKNOWN_NONZERO_PATTERN);
    MatMatMult(temp4, BST, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp5);
    finalize<fw>(temp5);
    MatMatMult(temp5, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp6);
    finalize<fw>(temp6);
    MatMatMult(temp6, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp7);
    finalize<fw>(temp7);
    MatMatMult(temp7, BS, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp8);
    finalize<fw>(temp8);
    MatAXPY(M, -1., temp8, UNKNOWN_NONZERO_PATTERN);

    destroy<fw>(LT);
    destroy<fw>(BST);
    destroy<fw>(temp1);
    destroy<fw>(temp2);
    destroy<fw>(temp3);
    destroy<fw>(temp4);
    destroy<fw>(temp5);
    destroy<fw>(temp6);
    destroy<fw>(temp7);
    destroy<fw>(temp8);

  }
}


void sbp_sat::x2::initialize_symbols(
  vv<std::size_t>       &F_symbols,
  vv<std::size_t>       &FT_symbols, 
  vv<std::size_t> const &interfaces,
  components      const &sbp) {

  constexpr std::size_t w = 1; constexpr std::size_t e = 2; 
  constexpr std::size_t s = 3; constexpr std::size_t n = 4; 

  for (std::size_t row = 0; row < sbp.n_blocks; ++row) {
    for (std::size_t col = 0; col < sbp.n_blocks; ++col) {
      std::size_t interface = interfaces[row][col];
      if (interface != 0 and row == col - 1) {  
        F_symbols[row][interface - 1] = n;
        F_symbols[col][interface - 1] = s;
        FT_symbols[interface - 1][row] = n;
        FT_symbols[interface - 1][col] = s;
      }
      else if (interface != 0 and row == col - 3) {  
        F_symbols[row][interface - 1] = e;
        F_symbols[col][interface - 1] = w;
        FT_symbols[interface - 1][row] = e;
        FT_symbols[interface - 1][col] = w;
      }
    }
  }
}

