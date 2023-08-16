#pragma once

#include <cmath>
#include <numbers>
#include <vector>
#include <tuple>
#include <iostream>
#include <array>
#include <optional>
#include <functional>
#include <iomanip>

#include "definitions.h"
#include "ranges.h"
#include "numerical.h"
#include "solve.h"

#include "linalg.h"
#include "components.h"
#include "connect.h"
#include "logging.h"
#include "local_solver.h"
#include "verify.h"

#include "compute_lambdaA_mkl.h"

#include "omp.h"

/* petsc headers */
#include "petscsys.h"
#include "petscmat.h" 
#include "petscviewer.h"
#include "petscvec.h"
#include "petscksp.h"

namespace sbp_sat { 

  using nat_t = std::size_t;
  using real_t = type::real_t;  
  
  template <typename t> 
  using vv = type::vv<t>;
  
  // using vector_t = std::vector<real_t>;
  using real_v = std::vector<real_t>;
  using size_v = std::vector<std::size_t>;
  using range_t = numerics::linrange<real_t>;

  using block_t = std::tuple<range_t, range_t>; 
  using block_v = std::vector<block_t>;

  using domain_t = std::tuple<real_t, real_t>;
  using domain_v = std::vector<domain_t>;

  using boundary_t = std::tuple<range_t, nat_t>;
  using boundary_tx2 = std::array<std::tuple<range_t, nat_t>, 4>;
  using boundary_vx2 = std::vector<boundary_tx2>;

  /*
  auto const π = std::numbers::pi_v<real_t>;
  auto static to_real_t = [](nat_t n){return static_cast<real_t>(n);};
  using petsc_matrix = Mat;
  */

  using petsc_vector = Vec;


namespace x2 {

  using ℤ = std::size_t;
  using ℝ = type::real_t;  

  // using namespace linalg;

  // auto constexpr fw = framework::petsc;

  template <typename T>
  using optref = std::optional<std::reference_wrapper<T>>;

  auto constexpr transposed = true;

  std::size_t constexpr dirichlet = 1;
  std::size_t constexpr neumann   = 2;
  std::size_t constexpr x         = 0;
  std::size_t constexpr y         = 1;
  std::size_t constexpr left      = 0;
  std::size_t constexpr right     = 1;

  /* */
  template <typename T>
  T analytical_solution(T &x, T &y, T const &cx=1., T const &cy=1.) {
      auto constexpr pi = std::numbers::pi_v<T>;
      return std::sin(cx * pi * x + cy * pi * y); 
  }

  template <typename T>
  T uxx(T &x, T &y, T const &cx=1., T const &cy=1.) {
      auto constexpr pi = std::numbers::pi_v<T>;
      return -(cx * cx) * (pi * pi) * std::sin(cx * pi * x + cy * pi * y); 
  }

  template <typename T>
  T uyy(T &x, T &y, T const &cx=1., T const &cy=1.) {
      auto constexpr pi = std::numbers::pi_v<T>;
      return -(cy * cy) * (pi * pi) * std::sin(cx * pi * x + cy * pi * y); 
  }

  template <typename T>
  T ux(T &x, T &y, T const &cx=1., T const &cy=1.) {
      auto constexpr pi = std::numbers::pi_v<T>;
      return cx * pi * std::cos(cx * pi * x + cy * pi * y); 
  }

  template <typename T>
  T uy(T &x, T &y, T const &cx=1., T const &cy=1.) {
      auto constexpr pi = std::numbers::pi_v<T>;
      return cy * pi * std::cos(cx * pi * x + cy * pi * y); 
  }

  void petsc_poisson(
    sbp_sat::real_v             &result,
    sbp_sat::domain_v     const &domain, 
    sbp_sat::block_t      const &block, 
    sbp_sat::boundary_tx2 const &boundaries); 

  void petsc_hybridized_poisson(std::size_t vl_n, std::size_t el_n);

  void write_m( 
    petsc_matrix       &m,
    block_t      const &block_x1, 
    block_t      const &block_x2, 
    std::array<std::size_t, 4> const bc);

  void write_d2_h1(
    petsc_matrix       &M, 
    real_v       const &h, 
    range_t      const &local,
    real_t       const spacing_square, 
    real_t       const coeff=1.); 

  void write_h1(
    petsc_matrix       &M, 
    real_v       const &h); 

  void write_boundary_x1_left(
    petsc_matrix       &M, 
    range_t      const &x1,
    range_t      const &x2,
    vector_t     const &bsx1,
    vector_t     const &hx2, 
    long double  const  β,
    long double  const  τ,
    int          const  order = 1);

  void write_boundary_x1_right(
    petsc_matrix       &M, 
    range_t      const &x1,
    range_t      const &x2,
    vector_t     const &bsx1,
    vector_t     const &hx2, 
    real_t       const  β,
    real_t       const  τ,
    int          const  order = 1); 

  void write_boundary_x2_left(
    petsc_matrix       &M, 
    range_t      const &x1,
    range_t      const &x2,
    vector_t     const &bsx2,
    vector_t     const &hx1, 
    real_t       const  β,
    real_t       const  τ,
    int          const  order = 1);

  void write_boundary_x2_right(
    petsc_matrix       &M, 
    range_t      const &x1,
    range_t      const &x2,
    vector_t     const &bsx2,
    vector_t     const &hx1, 
    real_t       const  β,
    real_t       const  τ,
    int          const  order = 1);

  void write_Ls(
    std::vector<petsc_matrix> &L,
    boundary_vx2 const &boundaries);

  void write_Lts(
    std::vector<sbp_sat::petsc_matrix> &L,
    std::vector<sbp_sat::petsc_matrix> &Lt); 

/* 

            (1, 0)           
            ______
           |      |
    (0, 0) |      | (0, 1)
           |______|

            (1, 1)

   _______      _______      _______      _______
  |\      |    |       |    |'      |    |,      |
  |  \    |    |       |    |  '    |    |  ,    |
  |       |    |    \  |    |    '  |    |    ,  |
  |_______|    |______\|    |______'|    |______,|
    
    (0,0)        (0,1)        (1,0)        (1,1)

*/

/* Set a partial Kronecker product i.e., lower or upper, for a single 
   diagonal matrix. Where the full Kronecker product would result in 
   entries on a single boundary row or column. */
template<nat_t major, nat_t minor> 
void boundary_diagonal(
  matrix<fw>       &dest,
  nat_t      const  size1,
  nat_t      const  size2,
  vector_t   const &h, 
  real_t     const  c = 1.) {

  nat_t row, col; 
  real_t val; 
  nat_t const size = size1 * size2;
  if constexpr (major == x and minor == left) {
    for (std::size_t i = 0; i < h.size(); ++i) {
      row = i;
      col = i;
      val = c * h[i];
      set_matrix_value<fw>(dest, row, col, val);
    }
  }
  else if constexpr (major == x and minor == right) {
    for (std::size_t i = 0; i < h.size(); ++i) {
      row = size - i - 1;
      col = size - i - 1;
      val = c * h[h.size() - i - 1];
      set_matrix_value<fw>(dest, row, col, val);
    } 
  }
  else if constexpr (major == y and minor == left) {
    for (std::size_t i = 0; i < h.size(); ++i) {
      row = i * size2;
      col = i * size2;
      val = c * h[i];
      set_matrix_value<fw>(dest, row, col, val);
    } 
  }
  else if constexpr (major == y and minor == right) {
    for (std::size_t i = 0; i < h.size(); ++i) { 
      col = size - (i * size2) - 1;
      val = c * h[h.size() - i - 1];
      row = size - (i * size2) - 1;
      set_matrix_value<fw>(dest, row, col, val);
    } 
  }
}


template<std::size_t major, std::size_t minor, bool transpose = false> 
void boundary_bs(
  matrix<fw>             & dest,
  nat_t            const   size1,
  nat_t            const   size2,
  vector_t         const & bs, 
  optref<vector_t  const>  h = {},
  real_t           const   c = 1.) {

  std::function<real_t(nat_t, nat_t)> entry;
  if (h) {
    vector_t tmp = *h;
    entry = [&bs, h = tmp, &c] (nat_t i, nat_t j) -> real_t {
      return c * h[i] * bs[j];};
  }
  else
    entry = [&bs, &c] (nat_t i, nat_t j) -> real_t {
        (void) i; return c * bs[j];};

  vector_t const & spacing = *h;
  nat_t row, col; 
  real_t val;
  nat_t const size = size1 * size2;
  if constexpr (major == x and minor == left and transpose) {
    for (std::size_t i = 0; i < spacing.size(); ++i) {
      col = i;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        row = i + (j * size1);
        set_matrix_value<fw>(dest, row, col, val);
      }
    }
  }
  else if constexpr (major == x and minor == left and not transpose) {
    for (std::size_t i = 0; i < size2; ++i) {
      row = i;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        col = i + (j * size1);
        set_matrix_value<fw>(dest, row, col, val);
      }
    }
  }
  else if constexpr (major == x and minor == right and transpose) {
    for (std::size_t i = 0; i < size2; ++i) {
      col = size - i - 1;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        row = size - i - (j * size1) - 1;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major == x and minor == right and not transpose) {
    for (std::size_t i = 0; i < size2; ++i) {
      row = size - i - 1;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        col = size - i - (j * size1) - 1;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major == y and minor == left and transpose) {
    for (std::size_t i = 0; i < size1; ++i) {
      col = i * size2;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        row = (i * size2) + j;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major == y and minor == left and not transpose) {
    for (std::size_t i = 0; i < size1; ++i) {
      row = i * size2;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        col = (i * size2) + j;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
  else if constexpr (major == y and minor == right and transpose) {
    for (std::size_t i = 0; i < size1; ++i) { 
      col = size - (i * size2) - 1;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        row = size - (i * size2) - j - 1;
        set_matrix_value<fw>(dest, row, col, val);
      }  
    } 
  }
  else if constexpr (major == y and minor == right and not transpose) {
    for (std::size_t i = 0; i < size1; ++i) { 
      row = size - (i * size2) - 1;
      for (std::size_t j = 0; j < bs.size(); ++j) {
        val = entry(i, j);
        col = size - (i * size2) - j - 1;
        set_matrix_value<fw>(dest, row, col, val);
      }
    } 
  }
}

template<std::size_t dim, std::size_t side, std::size_t order>
void add_boundary(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bs,
  vector_t     const &h, 
  real_t       const  β,
  real_t       const  τ) {

  matrix<fw> result;
  nat_t size = x1.size() * x2.size();

  if constexpr (order == neumann) {
    make_local_sparse_matrix<fw>(result, size, size, 3);
    boundary_diagonal<dim, side>(result, x1.size(), x2.size(), h, τ);
    boundary_bs<dim, side, transposed>(result, x1.size(), x2.size(), bs, h, -β);
  }
  else if constexpr (order == dirichlet) {
    matrix<fw> a;
    matrix<fw> b;
    make_local_sparse_matrix<fw>(a, size, size, 3);
    make_local_sparse_matrix<fw>(b, size, size, 3);
    boundary_diagonal<dim, side>(a, x1.size(), x2.size(), h);
    boundary_bs<dim, side, transposed>(a, x1.size(), x2.size(), bs, h, -1. / τ);
    boundary_bs<dim, side>(b, x1.size(), x2.size(), bs);
    finalize<fw>(a);
    finalize<fw>(b);
    matmul<fw>(a, b, result);
    destroy<fw>(a);
    destroy<fw>(b);
  }
  finalize<fw>(result);
  MatCompositeAddMat(M, result);
}

/*  Add boundary coefficients to a matrix operator for 2D volume points 
    given the dimension and side of the dimension that the boundary is 
    on for such a matrix. 

    Here, the matrix operator dimensions are organized by major and 
    minor axes, i.e., for major axis x, minor axis y, |x| = X, |y| = Y, 
    volume point (i, j) is denoted by row (i * X * Y) + (j * Y). 

    order selects either Neumann or Dirichlet boundaries
  
    Dirichlet:
    τ * H_x - β * H_x * BS_y' 
        
    Neumann:
    H_x * BS_y - (1 / τ) * H_x * BS_y * BS_y 

 */
template<std::size_t dim, std::size_t side>
void add_boundary(
  petsc_matrix       &M, 
  range_t      const &x1,
  range_t      const &x2,
  vector_t     const &bs,
  vector_t     const &h, 
  real_t       const  β,
  real_t       const  τ,
  std::size_t  const  order) {

  matrix<fw> result;
  nat_t size = x1.size() * x2.size();

  if (order == dirichlet) {
    make_local_sparse_matrix<fw>(result, size, size, 3);
    boundary_diagonal<dim, side>(result, x1.size(), x2.size(), h, τ);
    boundary_bs<dim, side, transposed>(result, x1.size(), x2.size(), bs, h, -β);
  }
  else if (order == neumann) {
    matrix<fw> a;
    matrix<fw> b;
    make_local_sparse_matrix<fw>(a, size, size, 3);
    make_local_sparse_matrix<fw>(b, size, size, 3);
    boundary_diagonal<dim, side>(a, x1.size(), x2.size(), h);
    boundary_bs<dim, side, transposed>(a, x1.size(), x2.size(), bs, h, -1. / τ);
    boundary_bs<dim, side>(b, x1.size(), x2.size(), bs);
    finalize<fw>(a);
    finalize<fw>(b);
    matmul<fw>(a, b, result);
    destroy<fw>(a);
    destroy<fw>(b);
  }
  finalize<fw>(result);
  MatCompositeAddMat(M, result);
}

void compute_solution(
  std::vector<petsc_vector>       &x,
  vv<petsc_solver>          const &A, 
  petsc_vector              const &b, 
  components                const &sbp);

void solve(
  KSP &A, 
  std::vector<petsc_vector> &b, 
  std::vector<petsc_vector> &x);

void make_F(
  components const &sbp, 
  std::vector<petsc_matrix> &F,
  vv<petsc_vector> &f,
  std::vector<std::vector<double>> &f_data);

void make_F_explicit(
  petsc_matrix           &F_explicit, 
  vv<petsc_vector> const &f,
  vv<std::size_t>  const &F_symbols,  
  components       const &sbp);

void make_F_sliced(
  std::vector<petsc_vector>       &F_sliced, 
  petsc_matrix              const &F_explicit,  
  components                const &sbp);

void make_F_textra(
  components const &sbp, 
  vv<petsc_matrix> &F);

void make_explicit_solvers(
  std::vector<KSP>                &solvers,
  std::vector<petsc_matrix> const &M,
  components const &sbp);

void fcompop(
  petsc_matrix &f, 
  petsc_matrix const &l, 
  petsc_matrix const &b,
  petsc_matrix const &h,
  real_t       const τ, 
  real_t       const β);

void make_scaled_vec(
  petsc_vector &v, 
  std::size_t   n,
  real_t        s = 1.);

void make_index_vec(std::vector<int> &v);

// Definitions located in trace_a.cpp. These include functions exclusive  
// to computing λA = D - FT * M \ F 

// Top-level functions
void compute_λA_localized(
  petsc_matrix                    &        λA,
  std::vector<KSP>          const &         M, 
  vv<petsc_vector>          const &         F,
  vv<std::size_t>           const &FT_symbols, 
  vv<std::size_t>           const & F_symbols, 
  vv<std::size_t>           const &interfaces, 
  components                const &       sbp);

// D Matrix functions
void compute_D(
  petsc_matrix          &D, 
  components      const &sbp, 
  vv<std::size_t> const &interfaces);

void make_interface_list(
  vv<std::size_t>       &interface_list, 
  vv<std::size_t>          const &F_symbols,
  vv<std::size_t>          const &FT_symbols,
  components               const &sbp); 


// M^-1 * F functions 
void compute_MF(
  vv<petsc_vector>       &x,
  vv<petsc_solver> const &m,
  vv<petsc_vector> const &f,
  components       const &sbp);

void initialize_MF(
  vv<petsc_vector>       &MF, 
  components       const &sbp);

void initialize_MF_reduced(
  vv<petsc_vector>       &MF, 
  components       const &sbp);

void destroy_MF(
  vv<petsc_vector> &MF);

void compute_MF_sliced(
  std::vector<petsc_vector>       &MF_sliced, 
  std::vector<KSP>          const &explicit_solvers, 
  std::vector<petsc_vector> const &F_sliced); 

void print_MF(
  vv<petsc_vector> const &MF, 
  vv<std::size_t> const &F_symbols, 
  components const &sbp);

// F^T * M^-1 * F functions
void compute_FTMF(
  vv<petsc_matrix>       &FTMF, 
  vv<petsc_vector> const &F, 
  vv<petsc_vector> const &MF, 
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp);

void ftmfcompop( // kernel used in compute_FTMF ie. (F^T) * (M^-1 F) 
  petsc_matrix                    &FTMF, 
  std::vector<petsc_vector> const &F, 
  std::vector<petsc_vector> const &MF);

void print_FTMF(
  vv<petsc_vector> const &F, 
  vv<petsc_vector> const &MF, 
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp);

void compute_λA( 
  petsc_matrix           &λA, 
  petsc_matrix     const &D, 
  vv<petsc_vector> const &F, 
  vv<petsc_vector> const &MF, 
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp); 

void compute_λA_reduced( 
  petsc_matrix           &λA, 
  petsc_matrix     const &D, 
  vv<petsc_vector> const &F, 
  vv<petsc_vector> const &MF, 
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp); 

void compute_λA_gemm( 
  petsc_matrix                    &λA, 
  petsc_matrix              const &D, 
  std::vector<petsc_matrix> const &F, 
  std::vector<petsc_matrix> const &MF, 
  vv<std::size_t>           const &F_symbols,
  vv<std::size_t>           const &FT_symbols,
  components                const &sbp); 

void compute_λA_gemm_flat( 
  petsc_matrix                    &λA, 
  petsc_matrix              const &D, 
  std::vector<petsc_matrix> const &MF, 
  vv<std::size_t>           const &indices, 
  components                const &sbp); 

void copy_MF(
  std::vector<petsc_matrix>       &MF_mat,
  vv<petsc_vector>          const &F_symbols);  

void initialize_λA(
  petsc_matrix     &λA, 
  components const &sbp);

// trace_b.cpp function 

void compute_sources(
    std::vector<petsc_matrix>             &F, 
    std::vector<range_t>            const &space,
    std::function<real_t(real_t, real_t)> const  f);

using boundary_functions = std::vector<
    std::function<real_t(real_t, real_t)>>;
using boundary_vectors = std::array<real_t, 4>;

void compute_boundary_solution(
    vv<petsc_vector>            &g, 
    std::vector<range_t>  const &grids,
    boundary_functions    const  bf, 
    boundary_vectors      const  b);

void compute_B(
    std::vector<petsc_matrix>       &  B, 
    components                const &sbp);

void compute_b1(
    petsc_matrix       & B, 
    petsc_matrix const & H, 
    real_t       const   τ, 
    petsc_matrix const & L, 
    real_t       const   β, 
    petsc_matrix const &BS,
    std::size_t  const   n);

void compute_b2(
    petsc_matrix       & B, 
    petsc_matrix const & H, 
    real_t       const   τ, 
    petsc_matrix const & L, 
    real_t       const   β, 
    petsc_matrix const &BS,
    std::size_t  const   n);

void compute_g(
    std::vector<petsc_vector>       &g, 
    std::vector<petsc_matrix> const &boundaries,
    vv<petsc_vector>          const &solutions,
    std::vector<petsc_matrix> const &sources, 
    vv<std::size_t>           const &boundary_type_map,
    vv<std::size_t>           const &boundary_data_map,
    components                const &sbp);

void compute_Mg(
  std::vector<petsc_vector>       &x,
  vv<petsc_solver>          const &M,
  std::vector<petsc_vector> const &g, 
  components                const &sbp);

void initialize_Mg(
  std::vector<petsc_vector>       &Mg, 
  components                const &sbp);

void destroy_Mg(
  std::vector<petsc_vector> &Mg);

void initialize_λb( 
  petsc_vector           &λb,
  components       const &sbp);

void compute_λb( 
  petsc_vector           &λb, 
  std::vector<petsc_matrix> const &F, 
  std::vector<petsc_vector> const &Mg, 
  vv<std::size_t>  const &F_symbols,
  components       const &sbp);


void make_M(
  petsc_matrix &M,
  components const &sbp, 
  std::array<std::size_t, 4> const &boundary);

void make_M_boundary(
  petsc_matrix &M,
  components const &sbp, 
  std::size_t const direction,
  std::size_t const boundary);

void initialize_symbols(
  vv<std::size_t> &F_Symbols,
  vv<std::size_t> &FT_Symbols, 
  vv<std::size_t> const &interfaces,
  components const &sbp);


// convergence.cpp
void copy_approx(
  petsc_vector &approx, 
  std::vector<petsc_vector> const &x, 
  components const &sbp); 

void compute_exact(
  petsc_vector &exact, 
  components const &sbp);

void print_convergence(
  petsc_vector &approx, 
  petsc_vector &exact, 
  components const &sbp); 

}; /* namespace sbp_sat::x2 */
}; /* namespace sbp_sat     */
