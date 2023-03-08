#include "hybrid_sbp_sat_2d.h"
/* Methods for computing the matrix portion of the trace variables. 
	
   g_δ - FT * M \ g 

   Requires known M, F, SBP-SAT components, and interfaces matrix, and 
   g (boundary terms.)
*/

// given a list of n ranges, generate n^2 matrices of source data: 
// (0, 0), (0, 1), ... (1, 0), (1, 1), ... (n-1, 0) ... (n-1, n-1)  
// assumed that f(x, y) = source_x(x, y) + source_y(x, y)
void sbp_sat::x2::compute_sources(
    std::vector<petsc_matrix>                   &F, 
    std::vector<range_t>                  const &grids,
    std::function<real_t(real_t, real_t)> const  f) {

    F.resize(grids.size() * grids.size());

    std::size_t i = 0;
    for (auto &x : grids) { for (auto &y : grids) {
        MatCreateSeqDense(PETSC_COMM_SELF, x.size(), y.size(), NULL, 
            &F[i]);
        for (auto xi = x.begin(); xi != x.end(); ++xi) { 
            for (auto yi = y.begin(); yi != y.end(); ++yi) {
                MatSetValue(F[i], xi.index, yi.index, f(*xi, *yi), 
                    INSERT_VALUES);
            }
        }
        finalize<fw>(F[i]);
        ++i;
    }}
    // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(F[8], PETSC_VIEWER_STDOUT_SELF);
}

// NOTE: Assumes number of grid points in x and y dims are equal. 
void sbp_sat::x2::compute_boundary_solution(
    vv<petsc_vector>            &g, 
    std::vector<range_t>  const &ranges,
    boundary_functions    const  bf, 
    boundary_vectors      const  b) {

    g.resize(4);
    for (auto &gg : g) gg.resize(ranges.size());

    for (std::size_t i = 0; i != 4; ++i) {
        for (std::size_t j = 0; j != ranges.size(); ++j) {
            auto range = ranges[j];            

            // std::cout << range.size() << std::endl;

            
            VecCreateSeq(PETSC_COMM_SELF, range.size(), &g[i][j]);

            auto inds = std::vector<int>(range.size());
            auto vals = std::vector<double>(range.size());
            for (auto e = range.begin(); e != range.end(); ++e) {
                inds[e.index] = e.index;
                vals[e.index] = (i < 2) 
                    ? bf[i](b[i], *e) : bf[i](*e, b[i]); 
            }

            VecSetValues(g[i][j], range.size(), &inds[0], &vals[0], 
                INSERT_VALUES);

            finalize<fw>(g[i][j]);
        }
    }
}

void sbp_sat::x2::compute_B(
    std::vector<petsc_matrix>       &  B, 
    components                const &sbp) {

    B.resize(8); // 4 for first-order 4 for second-order. Some may not be
                 // needed but it's easier to just compute them all for 
                 // now.

    compute_b1(B[0], sbp.hy, sbp.τ, sbp.lw, sbp.β, sbp.bsx, sbp.n);
    compute_b1(B[1], sbp.hy, sbp.τ, sbp.le, sbp.β, sbp.bsx, sbp.n);
    compute_b1(B[2], sbp.hx, sbp.τ, sbp.ls, sbp.β, sbp.bsy, sbp.n);
    compute_b1(B[3], sbp.hx, sbp.τ, sbp.ln, sbp.β, sbp.bsy, sbp.n);

    compute_b2(B[4], sbp.hy, sbp.τ, sbp.lw, sbp.β, sbp.bsx, sbp.n);
    compute_b2(B[5], sbp.hy, sbp.τ, sbp.le, sbp.β, sbp.bsx, sbp.n);
    compute_b2(B[6], sbp.hx, sbp.τ, sbp.ls, sbp.β, sbp.bsy, sbp.n);
    compute_b2(B[7], sbp.hx, sbp.τ, sbp.ln, sbp.β, sbp.bsy, sbp.n);
}

/* Compute H(a) * (τ * L(d)' - β * BS(a)' * L(d)') for particular 
   directional orientations of L and axes orientations of BS. For 2-D 
   EWNS this includes 

    East:  H_y * (τ * LE' - β * BS_x' * LE')
    West:  H_y * (τ * LW' - β * BS_x' * LW')
    South: H_x * (τ * LS' - β * BS_y' * LS')
    North: H_x * (τ * LN' - β * BS_y' * LN')  
*/
void sbp_sat::x2::compute_b1(
    petsc_matrix       & B, 
    petsc_matrix const & H, 
    real_t       const   τ, 
    petsc_matrix const & L, 
    real_t       const   β, 
    petsc_matrix const &BS,
    std::size_t  const   n) {

    petsc_matrix LT, BST, temp1, temp2, temp3;
    MatTranspose(L, MAT_INITIAL_MATRIX, &LT);
    MatTranspose(BS, MAT_INITIAL_MATRIX, &BST);
    finalize<fw>(LT);
    finalize<fw>(BST);

    // temp1 := β * BS_a^T
    MatCreateSeqAIJ(PETSC_COMM_SELF, n * n, n * n, n, nullptr, &temp1);
    finalize<fw>(temp1);    
    MatAXPY(temp1, β, BST, UNKNOWN_NONZERO_PATTERN);    

    // temp2 := temp1 * L_d^T 
    MatMatMult(temp1, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2);
    finalize<fw>(temp2);   

    // temp3 := τ * L_d^T 
    MatCreateSeqAIJ(PETSC_COMM_SELF, n * n, n, n, nullptr, &temp3);
    finalize<fw>(temp3);    
    MatAXPY(temp3, τ, LT, UNKNOWN_NONZERO_PATTERN);

    // temp3 := temp3 + (-1) * temp2 
    MatAXPY(temp3, -1., temp2, UNKNOWN_NONZERO_PATTERN);

    // B := H * temp3
    MatMatMult(H, temp3, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);
    finalize<fw>(B);

    destroy<fw>(LT);
    destroy<fw>(BST);
    destroy<fw>(temp1);
    destroy<fw>(temp2);
    destroy<fw>(temp3);
}


/* Compute H(a) * (β * L(d)' - 1/τ * BS(a')' * L(d)') for particular 
   directional of L and axes orientations of H and BS. For 2-D EWNS this 
   includes 

    South: H_x * (β * LS' - 1/τ * BS_y' * LS')
    North: H_x * (β * LN' - 1/τ * BS_y' * LN')  */
void sbp_sat::x2::compute_b2(
    petsc_matrix       & B, 
    petsc_matrix const & H, 
    real_t       const   τ, 
    petsc_matrix const & L, 
    real_t       const   β, 
    petsc_matrix const &BS,
    std::size_t  const   n) {

    petsc_matrix LT, BST, temp1, temp2, temp3;
    MatTranspose(L, MAT_INITIAL_MATRIX, &LT);
    MatTranspose(BS, MAT_INITIAL_MATRIX, &BST);
    finalize<fw>(LT);
    finalize<fw>(BST);

    // temp1 := 1. / τ * BS(a')^T
    MatCreateSeqAIJ(PETSC_COMM_SELF, n * n, n * n, n, nullptr, &temp1);
    finalize<fw>(temp1);    
    MatAXPY(temp1, 1. / τ, BST, UNKNOWN_NONZERO_PATTERN);    

    // temp2 := temp1 * L_d^T 
    MatMatMult(temp1, LT, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &temp2);
    finalize<fw>(temp2);   

    // temp3 := β * L_d^T 
    MatCreateSeqAIJ(PETSC_COMM_SELF, n * n, n, n, nullptr, &temp3);
    finalize<fw>(temp3);    
    MatAXPY(temp3, β, LT, UNKNOWN_NONZERO_PATTERN);

    // temp3 := temp3 + (-1) * temp2 
    MatAXPY(temp3, -1., temp2, UNKNOWN_NONZERO_PATTERN);

    // B := H * temp3
    MatMatMult(H, temp3, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);
    finalize<fw>(B);

    destroy<fw>(LT);
    destroy<fw>(BST);
    destroy<fw>(temp1);
    destroy<fw>(temp2);
    destroy<fw>(temp3);
}

void sbp_sat::x2::compute_g(
    std::vector<petsc_vector>       &g, 
    std::vector<petsc_matrix> const &boundaries,
    vv<petsc_vector>          const &solutions,
    std::vector<petsc_matrix> const &sources, 
    vv<std::size_t>           const &boundary_type_map,
    vv<std::size_t>           const &boundary_data_map,
    components                const &sbp) {

    // compute ith: (boundary * solution ...) - H_tilde sources[i]
    //                       bm       bm         stat     m2 block 
    //    v         v        v        v         v         v      
    // (b_LB_W * g_LB_W + b_LB_S * g_LB_S) - H_tilde * F_LB[:]

    petsc_vector temp2, temp3, temp5, temp6;
    auto v = std::vector<double>(sbp.n * sbp.n);
    auto mi = std::vector<int>(sbp.n * sbp.n);
    for (std::size_t i = 0; i != sbp.n * sbp.n; ++i) {
        mi[i] = i;
    } 

    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp2);
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp3);
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp5);
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &temp6);

    g.resize(sbp.n_blocks);

    for (std::size_t block = 0; block != sbp.n_blocks; ++block) {

        VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &g[block]);

        // 4 is known quantity as our blocks are rect. and orth.
        constexpr std::size_t faces = 4;
        for (std::size_t face = 0; face != faces; ++face) {

            auto b_type_index = boundary_type_map[block][face];
            if (b_type_index != 0) {
                
                auto ti = face + ((b_type_index - 1) * faces);
                auto boundary = boundaries[ti];

                auto di = boundary_data_map[block][face] - 1;
                auto solution = solutions[face][di];

                MatMult(boundary, solution, temp2);
                VecAXPY(g[block], 1., temp2); // g[block] += temp2 * 1.
            }
        }

        // Get the raw dense matrix data, effectively reshape to a vector
        MatGetValues(
            sources[block], sbp.n, &mi[0], sbp.n, &mi[0], &v[0]);
        VecSetValues(temp5, sbp.n * sbp.n, &mi[0], &v[0], INSERT_VALUES);
        finalize<fw>(temp5);
        MatMult(sbp.hl, temp5, temp6);
        VecAXPY(g[block], -1., temp6); // g[block] -= temp5 
        finalize<fw>(g[block]);
        // VecView(g[block], PETSC_VIEWER_STDOUT_SELF);
    }

    destroy<fw>(temp2);
    destroy<fw>(temp3);
    destroy<fw>(temp5);
    destroy<fw>(temp6);
}

void sbp_sat::x2::compute_Mg(
    std::vector<petsc_vector>       &Mg,
    std::vector<KSP>          const &M,
    std::vector<petsc_vector> const &g) {

    // #pragma omp parallel for
    for (std::size_t i = 0; i != M.size(); ++i) {
        KSPSolve(M[i], g[i], Mg[i]);
    }
}

void sbp_sat::x2::initialize_Mg(
  std::vector<petsc_vector>       &Mg, 
  components                const &sbp) {
    Mg.resize(sbp.n_blocks);
    for (auto &e : Mg) {
        VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n, &e);
    }
}

void sbp_sat::x2::destroy_Mg(
  std::vector<petsc_vector> &Mg) {
    for (std::size_t i = 0; i < Mg.size(); ++i) destroy<fw>(Mg[i]);
}


void sbp_sat::x2::initialize_λb( 
  petsc_vector           &λb,
  components       const &sbp) {
    VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n_interfaces, &λb);
}

void sbp_sat::x2::compute_λb( 
  petsc_vector           &λb, 
  std::vector<petsc_matrix> const &F, 
  std::vector<petsc_vector> const &Mg, 
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp) {

  std::size_t findex;
  auto v = std::vector<double>(sbp.n);

  auto read_indices = std::vector<int>(sbp.n);
  for (std::size_t i = 0; i != read_indices.size(); ++i) 
      read_indices[i] = i;
  auto write_indices = std::vector<int>(sbp.n);

  petsc_vector temp, n1; 
  VecCreateSeq(PETSC_COMM_SELF, sbp.n, &temp);
  VecCreateSeq(PETSC_COMM_SELF, sbp.n * sbp.n_interfaces, &n1); 

  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {

    // Update indices based on the super index
    for (std::size_t j = 0; j != write_indices.size(); ++j) 
      write_indices[j] = (i * sbp.n) + j;

    for (std::size_t j = 0; j != sbp.n_blocks; ++j) {

      findex = FT_symbols[i][j];
      if (findex != 0) {

        findex--;

        MatMult(F[findex], Mg[j], temp);
        VecGetValues(temp, sbp.n, &read_indices[0], &v[0]);
        VecSetValues(n1, sbp.n, &write_indices[0], &v[0], ADD_VALUES); 
      }
    }
  }

  // note: for some reason the block indices 0 and 2 on the solution are 
  //       way worse than the rest of the solution. unsure why. maybe 
  //       fix later. 
  
  VecAXPY(λb, -1., n1); // delta g is not used right now so just negate

  finalize<fw>(λb);

  destroy<fw>(temp);
  destroy<fw>(n1);
}

