#include "compute_lambda_a.h"

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include <iostream>

#define min(x,y) (((x) < (y)) ? (x) : (y))


auto compute_lambda_a(
  sparse_matrix_t              *lambdaA,
  sparse_matrix_t              *D,
  std::vector<sparse_matrix_t> &F,
  std::vector<real_t *> &MF,
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components                   &sbp) -> void {

  matrix_descr da;
  da.type = SPARSE_MATRIX_TYPE_GENERAL;

  double *ftmf, *FTMF;
  FTMF = (double *) mkl_malloc(sizeof(double) * sbp.n * sbp.n * F.size() * MF.size(), 64);
  sparse_status_t status;
  for (std::size_t i = 0; i != MF.size(); ++i) {
    for (std::size_t j = 0; j != F.size(); ++j) {
      
      ftmf = &FTMF[((i * MF.size()) + j) * sbp.n * sbp.n];
      status = mkl_sparse_d_mm(
        SPARSE_OPERATION_TRANSPOSE, 1., F[j], da, 
        SPARSE_LAYOUT_COLUMN_MAJOR, MF[i],
        sbp.n, sbp.n, 
        1., ftmf, sbp.n);
      mkl_sparse_status(status);
    }
  }

  std::size_t ilc = 0;
  std::vector<std::size_t> interface_list; 
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        if (F_symbols[k][i] > 0 && FT_symbols[j][k] > 0) {
          // findex = FT_symbols[i][k] - 1;
          // auto r = k % sbp.n_blocks_dim == 0 ? 0 
          //       : k % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          //       : 1;
          // mindex = (r * 4) + F_symbols[k][j] - 1;
          interface_list.push_back(i);
          interface_list.push_back(j);
          interface_list.push_back(k);
          ilc += 1;
          //std::cout << i << " " << j << std::endl;
        }
      }
    }
  }

  csr<double> fake{sbp.n_interfaces, sbp.n_interfaces};
  double *vals, *lat;
  vals = (double *) mkl_malloc(sizeof(double) * sbp.n * sbp.n * ilc, 64);

  std::size_t a, b, c, mindex, findex;
  for (std::size_t i = 0; i < interface_list.size(); i += 3) {
    a = interface_list[i];
    b = interface_list[i + 1];
    c = interface_list[i + 2];
    findex = FT_symbols[a][c] - 1;
    auto r = c % sbp.n_blocks_dim == 0 ? 0 
          : c % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          : 1;
    mindex = (r * 4) + F_symbols[c][b] - 1;
  
    ftmf = &FTMF[((mindex * MF.size()) + findex) * sbp.n * sbp.n];
    lat = &vals[i * sbp.n * sbp.n];
    for (std::size_t j = 0; j != sbp.n * sbp.n; ++j) {
      lat[j] -= ftmf[i]; // -= bc *-1 on D- FTMF
    }
  }

  //std::cout << "A" << std::endl;
  std::size_t lasta = -1, lastb = -1;
  for (std::size_t i = 0; i < interface_list.size(); i += 3) {
    a = interface_list[i];
    b = interface_list[i + 1];
    if (lasta != a || lastb != b) {
      fake(1., a, b); // generate bsr fake csr.
      //std::cout << a <<  " " <<  b << std::endl;
    }
    lasta = a; lastb =  b;
    
  }
  for(auto &e: fake.c) e += 1;
  for(auto &e: fake.r) e += 1;
  
  std::cout << ".....\n";

  sparse_matrix_t pbA;
  status = mkl_sparse_d_create_bsr(&pbA, 
    SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_COLUMN_MAJOR, 
    sbp.n_interfaces, sbp.n_interfaces, sbp.n, 
    &fake.r[0], &fake.r[1], &fake.c[0], vals);
  mkl_sparse_status(status);

  sparse_matrix_t dest; 
  status = mkl_sparse_convert_csr(pbA, 
    SPARSE_OPERATION_NON_TRANSPOSE, &dest);
  mkl_sparse_status(status);

  status =  mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, 
    *D, 1., dest, lambdaA);
  mkl_sparse_status(status);

  mkl_free(FTMF);
  mkl_free(vals); 
  mkl_sparse_destroy(pbA);
  mkl_sparse_destroy(dest);
  /*
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

          ftmf = &FTMF[((mindex * MF.size()) + findex) * sbp.n * sbp.n];

          memcpy( sbp.n * sbp.n);

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
  */
  return;
}

/* mkl_sparse_d_mm (
    const sparse_operation_t operation, 
    const double alpha, const sparse_matrix_t A, 
    const struct matrix_descr descr, 
    const sparse_layout_t layout, const double *B, 
    const MKL_INT columns, 
    const MKL_INT ldb,
     const double beta, double *C, 
     const MKL_INT ldc);*/


