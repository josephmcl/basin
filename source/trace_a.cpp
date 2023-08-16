#include "hybrid_sbp_sat_2d.h"
/* Methods for computing the matrix portion of the trace variables. 
	
   D - FT * M \ F 

   Requires known M, F, SBP-SAT components, and interfaces matrix.
*/

#include <iostream>
#include <chrono>

// Top level solver for D - FT * M \ F 
/*
void sbp_sat::x2::compute_λA_localized(
	petsc_matrix                    &        λA,
	std::vector<KSP>          const &         M, 
	vv<petsc_vector>          const &         F,
	vv<std::size_t>           const &FT_symbols, 
	vv<std::size_t>           const & F_symbols, 
 	vv<std::size_t>           const &interfaces, 
	components                const &       sbp) {

	auto MF = vv<petsc_vector>(4 * sbp.n_blocks, 
            std::vector<petsc_vector>(sbp.n)); 
  	initialize_MF(MF, sbp);
  	compute_MF(MF, M, F);

  	petsc_matrix D;
  	compute_D(D, sbp, interfaces);

  	compute_λA(λA, D, F, MF, F_symbols, FT_symbols, sbp);

  	std::cout << "done." << std::endl;

  	destroy_MF(MF);
  	destroy<fw>(D);
}
*/

void sbp_sat::x2::compute_D(
  petsc_matrix          &D, 
  components      const &sbp, 
  vv<std::size_t> const &interfaces) {

  /* D = [H1y * 2 * τ,                              ]
         [             H1y * 2 * τ,                 ]
         [                        , ...,            ]
         [                               H1y * 2 * τ] 

  NOTE: either H1x or H1y for corresponding interface 
        orientations (NS or EW). For now we use h1v because 
        h1x and h1y are identical in this current config. */

  std::size_t n_interfaces = 0;
  for (std::size_t row = 0; row != interfaces.size(); ++row) {
    for (std::size_t col = 0; col != interfaces.size(); ++col) {
      std::size_t interface = interfaces[row][col];
      if (interface != 0) n_interfaces += 1;
    }
  }

  std::size_t size = sbp.n * n_interfaces;
  MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 1, nullptr, &D);

  std::size_t interface, index;
  real_t value;
  for (std::size_t row = 0; row != interfaces.size(); ++row) {
    for (std::size_t col = 0; col != interfaces.size(); ++col) {
      interface = interfaces[row][col];
      if (interface != 0 and row == col - 1) {  // NS 
        index = (interface - 1) * sbp.n; 
        for (std::size_t k = 0; k != sbp.n; ++k) {
          value = sbp.h1v[k] * 2. * sbp.τ;
          MatSetValue(D, index + k,  index + k, value, ADD_VALUES);
        }
      }
      else if (interface != 0) {  // EW
        index = (interface - 1) * sbp.n; 
        for (std::size_t k = 0; k != sbp.n; ++k) {
          value = sbp.h1v[k] * 2. * sbp.τ;
          MatSetValue(D, index + k,  index + k, value, ADD_VALUES);
        }
      }
    }
  }
  finalize<fw>(D);
}

void sbp_sat::x2::initialize_MF(
  vv<petsc_vector>       &MF, 
  components       const &sbp) {

  for (std::size_t i = 0; i != 4 * sbp.n_blocks; ++i) {
    for (std::size_t j = 0; j != sbp.n; ++j) {
      VecCreateSeq(
        PETSC_COMM_SELF, 
        sbp.n * sbp.n, 
        &MF[i][j]);
    }
  } 
}

void sbp_sat::x2::initialize_MF_reduced(
  vv<petsc_vector>       &MF, 
  components       const &sbp) {

  for (std::size_t i = 0; i != 4 * sbp.n_blocks; ++i) {
    for (std::size_t j = 0; j != sbp.n; ++j) {
      VecCreateSeq(
        PETSC_COMM_SELF, 
        sbp.n * sbp.n, 
        &MF[i][j]);
    }
  } 
}

void sbp_sat::x2::destroy_MF(
  vv<petsc_vector> &MF) {

  for (std::size_t i = 0; i != MF.size(); ++i) {
    for (std::size_t j = 0; j != MF[i].size(); ++j) {
      destroy<fw>(MF[i][j]);
    }
  } 
}

void sbp_sat::x2::print_MF(
  vv<petsc_vector> const &MF, 
  vv<std::size_t> const &F_symbols, 
  components const &sbp) {

  int vn = sbp.n * sbp.n;
  auto vi = std::vector<int>(sbp.n * sbp.n);
  auto vy = std::vector<double>(sbp.n * sbp.n);
  for (auto i = 0; i < vn; ++i) vi[i] = i;

  petsc_matrix gMF;
  MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n_blocks * sbp.n * sbp.n,
    sbp.n_interfaces * sbp.n, sbp.n_interfaces * sbp.n, nullptr, 
    &gMF);

  std::cout << "M^-1xF size: " << MF.size() << " " 
    << MF[0].size() << std::endl;

  for (std::size_t i = 0; i < sbp.n_blocks; ++i) {
    for (std::size_t j = 0; j < sbp.n_interfaces; ++j) {
      if (F_symbols[i][j] > 0) {

        std::size_t m = i % sbp.n_blocks_dim == 0 ? 0 
          : i % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          : 1;
        std::size_t index = (m * 4) + F_symbols[i][j] - 1;
        std::size_t ii = i * sbp.n * sbp.n;
        std::size_t jj = j * sbp.n;
        
        std::cout << "retrieving M[" << i << "] and F[" 
          << F_symbols[i][j] << "] at index " << index << std::endl;

        for (std::size_t k = 0; k < sbp.n; ++k) {
          std::cout << k << std::endl; 
          VecGetValues(MF[index][k], vn, &vi[0], &vy[0]);
          std::cout << "got! " << std::endl; 
          for (std::size_t l = 0; l < sbp.n * sbp.n; ++l) {
            MatSetValue(gMF, ii + l, jj + k, vy[l], ADD_VALUES);
          }
        }
      }
    }
  }
  finalize<fw>(gMF);
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  MatView(gMF, PETSC_VIEWER_STDOUT_SELF);
}

// Given a vector of matrices M, and a 2-D vector of vectors F, solve 
// Mi x = Fi for all Mi in M and Fi in F. Store the result in X, aligned
// M-index major and F-index minor. 
void sbp_sat::x2::compute_MF(
  vv<petsc_vector>       &x,
  vv<petsc_solver> const &m,
  vv<petsc_vector> const &f,
  components       const &sbp) {
  
  std::size_t block_index, factor_index, slice_index, x_index, thread_index;

  std::size_t const limit = m.size() * f.size() * f[0].size();
  #pragma omp parallel for num_threads(sbp.n_threads) private(block_index, factor_index, slice_index, x_index, thread_index)
  for (std::size_t index = 0; index != limit; ++index) {
 
    block_index = index / (f.size() * f[0].size());
  	factor_index = (index % (f.size() * f[0].size())) / f[0].size();
    slice_index = (index % (f.size() * f[0].size())) - (factor_index * f[0].size()); 
    x_index = factor_index + (block_index * f.size());
    thread_index = static_cast<std::size_t>(omp_get_thread_num());

    KSPSolve(m[block_index][thread_index], f[factor_index][slice_index], x[x_index][slice_index]);
  }
}

void sbp_sat::x2::compute_MF_sliced(
  std::vector<petsc_vector>       &MF_sliced, 
  std::vector<KSP>          const &explicit_solvers, 
  std::vector<petsc_vector> const &F_sliced) {

  // #pragma omp parallel for num_threads(4) private(F_sliced)
  for (std::size_t i = 0; i != explicit_solvers.size(); ++i) {
    for (std::size_t j = 0; j != F_sliced.size(); ++j) { 
      if (j % explicit_solvers.size() == i) {
        KSPSolve(explicit_solvers[i], F_sliced[j], MF_sliced[j]);
      }
    }
  }
}

void sbp_sat::x2::compute_FTMF(
  vv<petsc_matrix>       &FTMF, 
  vv<petsc_vector> const &F, 
  vv<petsc_vector> const &MF, 
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp) {

  std::size_t findex;
  std::size_t mindex;

  // NOTE: j and k are both bound to block indices.
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {
          findex = FT_symbols[i][k] - 1;
          mindex = (k * 4) + F_symbols[k][j] - 1;
          // findex = 0;
          // mindex = 0;
          std::cout << "Compute FTMF super-index " << i 
           << ", " << j << " (" << FT_symbols[i][k] << "), (" 
           << F_symbols[k][j] << ")" << std::endl;
          ftmfcompop(FTMF[i][j], F[findex], MF[mindex]); 
        }
      }
      finalize<fw>(FTMF[i][j]);
    }
  }
}

void sbp_sat::x2::ftmfcompop(
  petsc_matrix                    &FTMF, 
  std::vector<petsc_vector> const &F, 
  std::vector<petsc_vector> const &MF) {

  double v;
  for (std::size_t i = 0; i != F.size(); ++i) {
    for (std::size_t j = 0; j != MF.size(); ++j) {
      v = 0.;
      VecTDot(F[i], MF[j], &v);
      MatSetValue(FTMF, i, j, v, ADD_VALUES); 
    }
  } 
}

void sbp_sat::x2::print_FTMF(
  vv<petsc_vector> const &F, 
  vv<petsc_vector> const &MF, 
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp) {

  petsc_matrix gFTMF;
  MatCreateSeqAIJ(PETSC_COMM_SELF, sbp.n_interfaces * sbp.n,
    sbp.n_interfaces * sbp.n, sbp.n_interfaces * sbp.n, nullptr, 
    &gFTMF);

  std::size_t findex;
  std::size_t mindex;
  // NOTE: j and k are both bound to block indices.
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {
          findex = FT_symbols[i][k] - 1;
          mindex = (k * 4) + F_symbols[k][j] - 1;
          // std::cout << "Compute GGGG FTMF super-index " << i 
          //  << ", " << j << " (" << FT_symbols[i][k] << "), (" 
          //  << F_symbols[k][j] << ")" << std::endl;

          double v;
          for (std::size_t ii = 0; ii != F[findex].size(); ++ii) {
            for (std::size_t jj = 0; jj != MF[mindex].size(); ++jj) {
              v = 0.;
              VecTDot(F[findex][ii], MF[mindex][jj], &v);
              MatSetValue(gFTMF, i * sbp.n + ii, j * sbp.n + jj, v, ADD_VALUES); 
            }
          } 

        }
      }
    }
  }

  finalize<fw>(gFTMF);
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  MatView(gFTMF, PETSC_VIEWER_STDOUT_SELF);
}

void sbp_sat::x2::compute_λA( 
  petsc_matrix           &λA, 
  petsc_matrix     const &D, 
  vv<petsc_vector> const &F, 
  vv<petsc_vector> const &MF, 
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp) {

  std::size_t findex;
  std::size_t mindex;
  // NOTE: j and k are both bound to block indices.

  #pragma omp parallel for collapse(2) private(findex, mindex) num_threads(sbp.n_threads)
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        // std::cout << omp_get_thread_num();
        if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {
          findex = FT_symbols[i][k] - 1;
          mindex = (k * 4) + F_symbols[k][j] - 1;
          // std::cout << "Compute GGGG FTMF super-index " << i 
          //  << ", " << j << " (" << FT_symbols[i][k] << "), (" 
          //  << F_symbols[k][j] << ")" << std::endl;

          double v;
          for (std::size_t ii = 0; ii != F[findex].size(); ++ii) {
            for (std::size_t jj = 0; jj != MF[mindex].size(); ++jj) {
              v = 0.;
              VecTDot(F[findex][ii], MF[mindex][jj], &v);
              MatSetValue(λA, i * sbp.n + ii, j * sbp.n + jj, v, ADD_VALUES); 
            }
          } 

        }
      }
    }
  }

  finalize<fw>(λA);

  MatAYPX(λA, -1, D, DIFFERENT_NONZERO_PATTERN);

  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  MatView(λA, PETSC_VIEWER_STDOUT_SELF);
}

void sbp_sat::x2::compute_λA_reduced( 
  petsc_matrix           &λA, 
  petsc_matrix     const &D, 
  vv<petsc_vector> const &F, 
  vv<petsc_vector> const &MF, 
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components       const &sbp) {

  std::size_t findex;
  std::size_t mindex;
  // NOTE: j and k are both bound to block indices.

  const std::size_t limit = sbp.n_interfaces * sbp.n_interfaces * sbp.n_blocks;
  std::size_t i, j, k;
  
  //std::size_t outer_counter = 0;
  //std::size_t inner_counter = 0;
  //std::size_t inner_nz_counter = 0;

  // std::cout << "limit :" << limit << std::endl;

  //int x = 0;
  //#pragma omp parallel for private(i, j, k, findex, mindex) num_threads(sbp.n_threads)
  //for (std::size_t index = 0; index != limit; ++index) {
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
  //  i = index / (sbp.n_interfaces * sbp.n_blocks);
  //  j = (index - (i * sbp.n_interfaces * sbp.n_blocks)) / sbp.n_blocks;
  //  k = (index - (j * sbp.n_blocks)) % sbp.n_blocks;
    
        if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {

          // outer_counter += 1;
          findex = FT_symbols[i][k] - 1;
          auto r = k % sbp.n_blocks_dim== 0 ? 0 
                : k % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
                : 1;
          mindex = (r * 4) + F_symbols[k][j] - 1;
          double v;

          for (std::size_t ii = 0; ii != F[findex].size(); ++ii) {
            for (std::size_t jj = 0; jj != MF[mindex].size(); ++jj) {
              v = 0.;
              VecTDot(F[findex][ii], MF[mindex][jj], &v);
              MatSetValue(λA, i * sbp.n + ii, j * sbp.n + jj, v, ADD_VALUES); 
            }
          } 
        }
      }
    }
  }

  finalize<fw>(λA);
  MatAYPX(λA, -1, D, DIFFERENT_NONZERO_PATTERN);
  // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  // MatView(λA, PETSC_VIEWER_STDOUT_SELF);
}

void sbp_sat::x2::compute_λA_gemm( 
  petsc_matrix                   &λA, 
  petsc_matrix              const &D, 
  std::vector<petsc_matrix> const &F, 
  std::vector<petsc_matrix> const &MF, 
  vv<std::size_t>           const &F_symbols,
  vv<std::size_t>           const &FT_symbols,
  components                const &sbp) {

  std::size_t findex;
  std::size_t mindex;
  // NOTE: j and k are both bound to block indices.

  const std::size_t limit = sbp.n_interfaces * sbp.n_interfaces * sbp.n_blocks;
  std::size_t i, j, k;
  
  //std::size_t outer_counter = 0;
  //std::size_t inner_counter = 0;
  //std::size_t inner_nz_counter = 0;

  // std::cout << "limit :" << limit << std::endl;

  //int x = 0;
  // #pragma omp parallel for private(i, j, k, findex, mindex) num_threads(sbp.n_threads)
  // #pragma omp parallel for private(i, j, k, findex, mindex) num_threads(sbp.n_threads)
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
    
        if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {

          // outer_counter += 1;
          findex = FT_symbols[i][k] - 1;
          auto r = k % sbp.n_blocks_dim == 0 ? 0 
                : k % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
                : 1;
          mindex = (r * 4) + F_symbols[k][j] - 1;

          std::vector<int> ind(sbp.n); 
          for (std::size_t i = 0; i != sbp.n; ++i) ind[i] = i;

          petsc_matrix res;
          MatMatMult(F[findex], MF[mindex], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &res);
          finalize<fw>(res);

          // int a,b,c,d;
          // MatGetSize(res, &a, &b);
          // std::cout << a << "-- " << b << std::endl;

          // std::cout << i << " " << j << std::endl;

          std::vector<int> allm(sbp.n); make_index_vec(allm);
          std::vector<int> alln(sbp.n) ; make_index_vec(alln);
          std::vector<double> data(sbp.n * sbp.n);
          MatGetValues(res, sbp.n, &allm[0], sbp.n, &alln[0], &data[0]);
          for (std::size_t b = 0; b != alln.size(); ++b) {
            allm[b] += i * sbp.n;
            alln[b] += j * sbp.n;
          }
          MatSetValues(λA, sbp.n, &allm[0], sbp.n, &alln[0], &data[0], ADD_VALUES);
        }
      }
    }
  }

  finalize<fw>(λA);
  MatAYPX(λA, -1, D, DIFFERENT_NONZERO_PATTERN);
  // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  // MatView(λA, PETSC_VIEWER_STDOUT_SELF);
}


void sbp_sat::x2::make_interface_list(
  vv<std::size_t>       &interface_list, 
  vv<std::size_t>          const &F_symbols,
  vv<std::size_t>          const &FT_symbols,
  components               const &sbp) {  

  std::size_t findex;
  std::size_t mindex;

  interface_list.resize(sbp.n_interfaces);
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        if (F_symbols[k][i] > 0 && FT_symbols[j][k] > 0) {
          interface_list[i].push_back(j);

          findex = FT_symbols[i][k] - 1;
          auto r = k % sbp.n_blocks_dim == 0 ? 0 
                : k % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
                : 1;
          mindex = (r * 4) + F_symbols[k][j] - 1;

          std::cout << findex << " " << mindex << std::endl;

        }
      }
    }
  }
}

void sbp_sat::x2::compute_λA_gemm_flat( 
  petsc_matrix                    &λA, 
  petsc_matrix              const &D, 
  std::vector<petsc_matrix> const &MF, 
  vv<std::size_t>           const &indices, 
  components                const &sbp) {

  std::size_t findex;
  std::size_t mindex;
  // NOTE: j and k are both bound to block indices.

  const std::size_t limit = sbp.n_interfaces * sbp.n_interfaces * sbp.n_blocks;
  std::size_t i, j, k;
  
  //std::size_t outer_counter = 0;
  //std::size_t inner_counter = 0;
  //std::size_t inner_nz_counter = 0;

  // std::cout << "limit :" << limit << std::endl;

  //int x = 0;
  // #pragma omp parallel for private(i, j, k, findex, mindex) num_threads(sbp.n_threads)
  // #pragma omp parallel for private(i, j, k, findex, mindex) num_threads(sbp.n_threads)

  for (std::size_t i = 0; i != indices.size(); ++i) {
    
    // findex[]

  }

  finalize<fw>(λA);
  MatAYPX(λA, -1, D, DIFFERENT_NONZERO_PATTERN);
  // PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  // MatView(λA, PETSC_VIEWER_STDOUT_SELF);

}

void 
sbp_sat::x2::copy_MF(
  std::vector<petsc_matrix>       &MF_mat,
  vv<petsc_vector>          const &MF_vec) {

  MF_mat.resize(MF_vec.size());

  int m, n, col;
  VecGetSize(MF_vec[0][0], &m);
  std::vector<int> rows(m);
  for (std::size_t i = 0; i != static_cast<std::size_t>(m); ++i) {
    rows[i] = i;
  }

  n = MF_vec[0].size();

  double *data;
  for (std::size_t i = 0; i != MF_vec.size(); ++i) {

    // make petsc matrix
    MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, n, nullptr, &MF_mat[i]);         

    // over each column 
    for (std::size_t j = 0; j != static_cast<std::size_t>(n); ++j) {

      col = j;

      // copy column vector data
      VecGetArray(MF_vec[i][j], &data); 

      // write into new matrix
      MatSetValues(MF_mat[i], m, &rows[0], 1, &col, data, INSERT_VALUES);
    }

    finalize<fw>(MF_mat[i]);
  }
} 

void sbp_sat::x2::initialize_λA(
  petsc_matrix     &λA, 
  components const &sbp) {

  MatCreateSeqAIJ(
    PETSC_COMM_SELF, 
    sbp.n_interfaces * sbp.n,
    sbp.n_interfaces * sbp.n, 
    sbp.n * (sbp.n_interfaces / 2), 
    nullptr, 
    &λA);
}

