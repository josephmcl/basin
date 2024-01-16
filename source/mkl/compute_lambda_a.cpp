#include "compute_lambda_a.h"

#include <iostream>
#include <cmath>

auto compute_lambda_a(
  double              *lambdaA,
  sparse_matrix_t              *D,
  std::vector<sparse_matrix_t> &F,
  std::vector<real_t *> &MF,
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components                   &sbp) -> void {

  
  matrix_descr da;
  da.type = SPARSE_MATRIX_TYPE_GENERAL;

  // std::cout << "sz " <<  F.size() * MF.size() << std::endl;
  double *ftmf, *FTMF, *lat;
  std::size_t sz = sizeof(double) * sbp.n * sbp.n * F.size() * MF.size();
  FTMF = (double *) mkl_malloc(sz, 64);
  std::memset(FTMF, 0, sz);
  sparse_status_t status;

  #pragma omp parallel for private(ftmf) collapse(2) num_threads(sbp.n_threads)
  for (std::size_t i = 0; i != MF.size(); ++i) {
    for (std::size_t j = 0; j != F.size(); ++j) {
      ftmf = &FTMF[(i * F.size() + j) * sbp.n * sbp.n];
      status = mkl_sparse_d_mm(
        SPARSE_OPERATION_NON_TRANSPOSE, 1., F[j], da, 
        SPARSE_LAYOUT_COLUMN_MAJOR, MF[i],
        sbp.n, sbp.n * sbp.n, 
        1., ftmf, sbp.n);
      mkl_sparse_status(status);
    }
  }


  std::size_t ilc = 0;
  std::vector<std::size_t> interface_list; 
  std::vector<std::size_t> iflk; 
  bool dirty = false;
  std::size_t slu = 0;
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      dirty = false;
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {
          if (!dirty) {
            slu += 1;
            dirty = true;
            iflk.push_back(i);
            iflk.push_back(j);
          }
          // findex = FT_symbols[i][k] - 1;
          // auto r = k % sbp.n_blocks_dim == 0 ? 0 
          //       : k % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          //       : 1;
          // mindex = (r * 4) + F_symbols[k][j] - 1;
          interface_list.push_back(i);
          interface_list.push_back(j);
          interface_list.push_back(k);
          interface_list.push_back(slu - 1);

          ilc += 1;
          //std::cout << i << " " << j << std::endl;
          
        }
      }
    }
  }

  // Copy each block into a new array, sometimes wit
  std::size_t a, b, c, d, mindex, findex;
  #pragma omp parallel for private(a, b, c, d, mindex, findex, ftmf, lat) num_threads(sbp.n_threads)
  for (std::size_t i = 0; i < interface_list.size(); i += 4) {
    a = interface_list[i];
    b = interface_list[i + 1];
    c = interface_list[i + 2];
    d = interface_list[i + 3];
    findex = FT_symbols[a][c] - 1;
    auto r = c % sbp.n_blocks_dim == 0 ? 0 
          : c % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          : 1;
    mindex = (r * 4) + F_symbols[c][b] - 1;
    // std::cout << i / 4 << " " << mindex << " " << findex << 
    //   " " << mindex * F.size() + findex <<std::endl;
    ftmf = &FTMF[(mindex * F.size() + findex) * sbp.n * sbp.n];
    
    std::size_t kl;
    
    for (std::size_t j = 0; j != sbp.n; ++j) {
      lat = &lambdaA[((j + (a * sbp.n)) * sbp.n_interfaces * sbp.n) + (b * sbp.n)];
      for (std::size_t k = 0; k != sbp.n; ++k) {
        lat[k] -= ftmf[j * sbp.n + k]; // -= bc *-1 on D- FTMF
        // kl = ((j + (a * sbp.n)) * sbp.n_interfaces * sbp.n) + (b * sbp.n) + k;
      }
    }
    // std::cout << std::endl;
  }

   /*
  std::cout << slu << std::endl;
  for (std::size_t s = 0; s != F.size() * MF.size(); ++s) {
    std::cout << s << ": ";
    for (std::size_t i = 0; i != sbp.n * sbp.n; ++i) {
      std::cout << FTMF[s * sbp.n * sbp.n + i] << " ";
    }
    std::cout <<  std::endl;
    std::cout <<  std::endl;
  }
 
  std::cout <<  std::endl;
  std::cout <<  std::endl;

    */

  /*
  for (std::size_t i = 0; i !=  sbp.n_interfaces * sbp.n; ++i) {
    for (std::size_t j = 0; j !=  sbp.n_interfaces * sbp.n; ++j) {
      std::cout << lambdaA[(i * sbp.n_interfaces * sbp.n) + j] << " ";
      //if (std::isnan(lambdaA[(i * sbp.n_interfaces * sbp.n) + j])) {
      //  std::cout << "!";
      //}
    } 
    std::cout <<  std::endl;
  } 

  std::cout << "??" << std::endl;
  */


  std::size_t v;
  #pragma omp parallel for private(v) num_threads(sbp.n_threads)
  for (std::size_t j = 0; j != sbp.n_interfaces * sbp.n; ++j) {
    v = j % sbp.n;
    lambdaA[(j * sbp.n_interfaces * sbp.n) + j] += sbp.h1v[v] * 2. * sbp.τ;
    // lambdaA[(j * sbp.n_interfaces * sbp.n) + j] += 1;
  }
  

  mkl_free(FTMF);

  /*


  sparse_matrix_t pbA;
  status = mkl_sparse_d_create_coo(&pbA, 
    SPARSE_INDEX_BASE_ZERO, 
    sbp.n * sbp.n_interfaces, sbp.n * sbp.n_interfaces, 
    sbp.n * sbp.n * slu, 
    &rowz[0], &colz[0], &vals[0]);
  mkl_sparse_status(status);
  
  //std::cout << fake.c.size() <<  " = " << sbp.n_interfaces << std::endl;
  //for (auto &e: fake.c)
  //  std::cout << e << " ";

  std::cout << "???" << std::endl;

  sparse_matrix_t dest; 
  status = mkl_sparse_convert_csr(pbA, 
    SPARSE_OPERATION_NON_TRANSPOSE, &dest);
  mkl_sparse_status(status);

  mkl_free(FTMF);
  mkl_free(vals); 
  mkl_free(rowz); 
  mkl_free(colz); 

  std::cout << "????" << std::endl;

  // mkl_sparse_order(dest);

  // std::cout << "?" << std::endl;

  
  sparse_index_base_t indexing;
    MKL_INT rows, cols, *rowst, *rowe, *coli, *ia;
    real_t *vz;
    status = mkl_sparse_d_export_csr(
      dest, &indexing, &rows, &cols, &rowst, &rowe, &coli, &vz);
    mkl_sparse_status(status);
    
    ia = (MKL_INT *) MKL_malloc(sizeof(MKL_INT) * ((sbp.n * sbp.n_interfaces) + 1), 64);
    std::memcpy(&ia[0], &rowst[0], sizeof(MKL_INT) * sbp.n * sbp.n_interfaces);
    ia[sbp.n * sbp.n_interfaces] = rowe[sbp.n * sbp.n_interfaces - 1];

    std::cout << rows << std::endl;
    std::cout << cols << std::endl;

    status = mkl_sparse_d_export_csr(
      *D, &indexing, &rows, &cols, &rowst, &rowe, &coli, &vz);
    mkl_sparse_status(status);

    std::cout << rows << std::endl;
    std::cout << cols << std::endl;

  
  status =  mkl_sparse_d_add(
    SPARSE_OPERATION_NON_TRANSPOSE, 
    *D, 1., dest, lambdaA);
  mkl_sparse_status(status);
  

  // std::cout << "?????" << std::endl;

  // status = mkl_sparse_qr_reorder(*lambdaA, da);
  // mkl_sparse_status(status);

  std::cout << "?????" << std::endl;
  
  mkl_sparse_destroy(pbA);
  mkl_sparse_destroy(dest);
  */
  return;
}

/*
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

  // std::cout << "sz " <<  F.size() * MF.size() << std::endl;
  double *ftmf, *FTMF;
  std::size_t sz = sizeof(double) * sbp.n * sbp.n * F.size() * MF.size();
  FTMF = (double *) mkl_malloc(sz, 64);
  std::memset(FTMF, 0, sz);
  sparse_status_t status;
  for (std::size_t i = 0; i != MF.size(); ++i) {
    for (std::size_t j = 0; j != F.size(); ++j) {

      ftmf = &FTMF[(i * F.size() + j) * sbp.n * sbp.n];
      status = mkl_sparse_d_mm(
        SPARSE_OPERATION_NON_TRANSPOSE, 1., F[j], da, 
        SPARSE_LAYOUT_COLUMN_MAJOR, MF[i],
        sbp.n, sbp.n * sbp.n, 
        1., ftmf, sbp.n);
      mkl_sparse_status(status);
    }
  }


  std::size_t ilc = 0;
  std::vector<std::size_t> interface_list; 
  std::vector<std::size_t> iflk; 
  bool dirty = false;
  std::size_t slu = 0;
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      dirty = false;
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {
          if (!dirty) {
            slu += 1;
            dirty = true;
            iflk.push_back(i);
            iflk.push_back(j);
          }
          // findex = FT_symbols[i][k] - 1;
          // auto r = k % sbp.n_blocks_dim == 0 ? 0 
          //       : k % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          //       : 1;
          // mindex = (r * 4) + F_symbols[k][j] - 1;
          interface_list.push_back(i);
          interface_list.push_back(j);
          interface_list.push_back(k);
          interface_list.push_back(slu - 1);

          ilc += 1;
          //std::cout << i << " " << j << std::endl;
          
        }
      }
    }
  }

  // std::cout << "??" << std::endl;
  
  double *vals, *lat; 
  MKL_INT *rowz, *colz, *rz, *cz;
  sz = sizeof(double) * sbp.n * sbp.n * slu;
  vals = (double *) mkl_malloc(sz, 64);
  rowz = (MKL_INT *) mkl_malloc(sz, 64);
  colz = (MKL_INT *) mkl_malloc(sz, 64);
  std::memset(vals, 0, sz);
  std::memset(rowz, 0, sz);
  std::memset(colz, 0, sz);

  // std::cout << "??" << std::endl;

  // Copy each block into a new array, sometimes wit
  std::size_t a, b, c, d, mindex, findex;
  for (std::size_t i = 0; i < interface_list.size(); i += 4) {
    a = interface_list[i];
    b = interface_list[i + 1];
    c = interface_list[i + 2];
    d = interface_list[i + 3];
    findex = FT_symbols[a][c] - 1;
    auto r = c % sbp.n_blocks_dim == 0 ? 0 
          : c % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          : 1;
    mindex = (r * 4) + F_symbols[c][b] - 1;
    // std::cout << i / 4 << " " << mindex << " " << findex << 
    //   " " << mindex * F.size() + findex <<std::endl;
    ftmf = &FTMF[(mindex * F.size() + findex) * sbp.n * sbp.n];
    lat = &vals[d * sbp.n * sbp.n];
    for (std::size_t j = 0; j != sbp.n * sbp.n; ++j) {
      lat[j] -= ftmf[j]; // -= bc *-1 on D- FTMF
    }
  }

  for (std::size_t i = 0; i < iflk.size(); i += 2) {
    a = iflk[i];
    b = iflk[i + 1];
    cz = &colz[(i/2) * sbp.n * sbp.n];
    rz = &rowz[(i/2) * sbp.n * sbp.n];
    for (std::size_t j = 0; j != sbp.n; ++j) {
      for (std::size_t k = 0; k != sbp.n; ++k) {
        cz[j * sbp.n + k] = a * sbp.n + j;
        rz[j * sbp.n + k] = b * sbp.n + k;
      }
    }    
  }

  std::cout << "??" << std::endl;

  // for(auto &e: fake.c) e += 1;
  // for(auto &e: fake.r) e += 1;
  
  sparse_matrix_t pbA;
  status = mkl_sparse_d_create_coo(&pbA, 
    SPARSE_INDEX_BASE_ZERO, 
    sbp.n * sbp.n_interfaces, sbp.n * sbp.n_interfaces, 
    sbp.n * sbp.n * slu, 
    &rowz[0], &colz[0], &vals[0]);
  mkl_sparse_status(status);
  
  //std::cout << fake.c.size() <<  " = " << sbp.n_interfaces << std::endl;
  //for (auto &e: fake.c)
  //  std::cout << e << " ";

  std::cout << "???" << std::endl;

  sparse_matrix_t dest; 
  status = mkl_sparse_convert_csr(pbA, 
    SPARSE_OPERATION_NON_TRANSPOSE, &dest);
  mkl_sparse_status(status);

  mkl_free(FTMF);
  mkl_free(vals); 
  mkl_free(rowz); 
  mkl_free(colz); 

  std::cout << "????" << std::endl;

  // mkl_sparse_order(dest);

  // std::cout << "?" << std::endl;

  
  sparse_index_base_t indexing;
    MKL_INT rows, cols, *rowst, *rowe, *coli, *ia;
    real_t *vz;
    status = mkl_sparse_d_export_csr(
      dest, &indexing, &rows, &cols, &rowst, &rowe, &coli, &vz);
    mkl_sparse_status(status);
    
    ia = (MKL_INT *) MKL_malloc(sizeof(MKL_INT) * ((sbp.n * sbp.n_interfaces) + 1), 64);
    std::memcpy(&ia[0], &rowst[0], sizeof(MKL_INT) * sbp.n * sbp.n_interfaces);
    ia[sbp.n * sbp.n_interfaces] = rowe[sbp.n * sbp.n_interfaces - 1];

    std::cout << rows << std::endl;
    std::cout << cols << std::endl;

    status = mkl_sparse_d_export_csr(
      *D, &indexing, &rows, &cols, &rowst, &rowe, &coli, &vz);
    mkl_sparse_status(status);

    std::cout << rows << std::endl;
    std::cout << cols << std::endl;

  
  status =  mkl_sparse_d_add(
    SPARSE_OPERATION_NON_TRANSPOSE, 
    *D, 1., dest, lambdaA);
  mkl_sparse_status(status);
  

  // std::cout << "?????" << std::endl;

  // status = mkl_sparse_qr_reorder(*lambdaA, da);
  // mkl_sparse_status(status);

  std::cout << "?????" << std::endl;
  
  mkl_sparse_destroy(pbA);
  mkl_sparse_destroy(dest);
  
  return;
}
*/

/*
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

  // std::cout << "sz " <<  F.size() * MF.size() << std::endl;
  double *ftmf, *FTMF;
  std::size_t sz = sizeof(double) * sbp.n * sbp.n * F.size() * MF.size();
  FTMF = (double *) mkl_malloc(sz, 64);
  std::memset(FTMF, 0, sz);
  sparse_status_t status;
  for (std::size_t i = 0; i != MF.size(); ++i) {
    for (std::size_t j = 0; j != F.size(); ++j) {

      // std::cout << i * F.size() + j << std::endl;
      // TODO HERE
      ftmf = &FTMF[(i * F.size() + j) * sbp.n * sbp.n];
      status = mkl_sparse_d_mm(
        SPARSE_OPERATION_NON_TRANSPOSE, 1., F[j], da, 
        SPARSE_LAYOUT_COLUMN_MAJOR, MF[i],
        sbp.n, sbp.n * sbp.n, 
        1., ftmf, sbp.n);
      mkl_sparse_status(status);
    }
  }

  // std::cout << "??" << std::endl;

  std::size_t ilc = 0;
  std::vector<std::size_t> interface_list; 
  std::vector<std::size_t> iflk; 
  bool dirty = false;
  std::size_t slu = 0;
  for (std::size_t i = 0; i != sbp.n_interfaces; ++i) {
    for (std::size_t j = 0; j != sbp.n_interfaces; ++j) {
      dirty = false;
      for (std::size_t k = 0; k != sbp.n_blocks; ++k) {
        if (FT_symbols[i][k] > 0 && F_symbols[k][j] > 0) {
          if (!dirty) {
            slu += 1;
            dirty = true;
            iflk.push_back(i);
            iflk.push_back(j);
          }
          // findex = FT_symbols[i][k] - 1;
          // auto r = k % sbp.n_blocks_dim == 0 ? 0 
          //       : k % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          //       : 1;
          // mindex = (r * 4) + F_symbols[k][j] - 1;
          interface_list.push_back(i);
          interface_list.push_back(j);
          interface_list.push_back(k);
          interface_list.push_back(slu - 1);

          ilc += 1;
          //std::cout << i << " " << j << std::endl;
          
        }
      }
    }
  }

  // std::cout << "??" << std::endl;
  
  csr<double> fake{sbp.n_interfaces, sbp.n_interfaces};
  double *vals, *lat;
  sz = sizeof(double) * sbp.n * sbp.n * (slu + 2);
  vals = (double *) mkl_malloc(sz, 64);
  std::memset(vals, 0, sz);

  // std::cout << "??" << std::endl;

  // Copy each block into a new array, sometimes wit
  std::size_t a, b, c, d, mindex, findex;
  for (std::size_t i = 0; i < interface_list.size(); i += 4) {
    a = interface_list[i];
    b = interface_list[i + 1];
    c = interface_list[i + 2];
    d = interface_list[i + 3];
    findex = FT_symbols[a][c] - 1;
    auto r = c % sbp.n_blocks_dim == 0 ? 0 
          : c % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          : 1;
    mindex = (r * 4) + F_symbols[c][b] - 1;
    // std::cout << i / 4 << " " << mindex << " " << findex << 
    //   " " << mindex * F.size() + findex <<std::endl;
    ftmf = &FTMF[(mindex * F.size() + findex) * sbp.n * sbp.n];
    lat = &vals[d * sbp.n * sbp.n];
    for (std::size_t j = 0; j != sbp.n * sbp.n; ++j) {
      lat[j] -= ftmf[j]; // -= bc *-1 on D- FTMF
    }
  }

  // std::cout << "??" << std::endl;

  //std::cout << "A" << std::endl;
  std::size_t lasta = -1, lastb = -1;
  for (std::size_t i = 0; i < iflk.size(); i += 2) {
    a = iflk[i];
    b = iflk[i + 1];
    if (lasta != a || lastb != b) {
      // std::cout << a <<  " " <<  b << std::endl;
      fake(1., a, b); // generate bsr fake csr.
      
    }
    lasta = a; lastb =  b;
    
  }
  
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

  status =  mkl_sparse_d_add(
    SPARSE_OPERATION_NON_TRANSPOSE, 
    *D, 1., dest, lambdaA);
  mkl_sparse_status(status);

  mkl_free(FTMF);
  mkl_free(vals); 
  mkl_sparse_destroy(pbA);
  mkl_sparse_destroy(dest);
  
  return;
}
*/

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

/* mkl_sparse_d_mm (
    const sparse_operation_t operation, 
    const double alpha, const sparse_matrix_t A, 
    const struct matrix_descr descr, 
    const sparse_layout_t layout, const double *B, 
    const MKL_INT columns, 
    const MKL_INT ldb,
     const double beta, double *C, 
     const MKL_INT ldc);*/


