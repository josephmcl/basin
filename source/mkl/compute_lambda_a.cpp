#include "compute_lambda_a.h"

#include <iostream>
#include <cmath>

// #include <papi.h>

auto compute_lambda_a_mpi(
  double              *lambdaA,
  sparse_matrix_t              *D,
  std::vector<sparse_matrix_t> &F,
  std::vector<real_t *> &MF,
  vv<std::size_t>  const &F_symbols,
  vv<std::size_t>  const &FT_symbols,
  components                   &sbp) -> void {

  
  matrix_descr da;
  da.type = SPARSE_MATRIX_TYPE_GENERAL;

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
        }
      }
    }
  }
  

  // Copy each block into a new array, sometimes wit
  std::size_t a, b, c, d, mindex, findex;
  //#pragma omp parallel for private(a, b, c, d, mindex, findex, ftmf, lat) // num_threads(sbp.n_threads)
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
    
    // std::size_t kl;
    
    

    for (std::size_t j = 0; j != sbp.n; ++j) {
      lat = &lambdaA[((j + (a * sbp.n)) * sbp.n_interfaces * sbp.n) + (b * sbp.n)];
      for (std::size_t k = 0; k != sbp.n; ++k) {
        lat[k] -= ftmf[j * sbp.n + k]; // -= bc *-1 on D- FTMF
        // kl = ((j + (a * sbp.n)) * sbp.n_interfaces * sbp.n) + (b * sbp.n) + k;
      }
    }
  }

  std::size_t v;
  // #pragma omp parallel for private(v) num_threads(sbp.n_threads)
  for (std::size_t j = 0; j != sbp.n_interfaces * sbp.n; ++j) {
    v = j % sbp.n;
    lambdaA[(j * sbp.n_interfaces * sbp.n) + j] += sbp.h1v[v] * 2. * sbp.τ;
    // lambdaA[(j * sbp.n_interfaces * sbp.n) + j] += 1;
  }
  

  mkl_free(FTMF);
  return;
}


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
  for (std::size_t i = 0; i < sbp.n * sbp.n * F.size() * MF.size(); ++i) {
    FTMF[i] = 0;
  }
  sparse_status_t status;


  // data layout: 
  // [ ft block -- 4                ]
  // [ m block -- 3         ]
  // [f block -- 4]
  // #pragma omp parallel for private(ftmf) collapse(2) num_threads(sbp.n_threads)
  for (std::size_t i = 0; i != F.size(); ++i) {
    for (std::size_t j = 0; j != MF.size(); ++j) {
      ftmf = &FTMF[(i * MF.size() + j) * sbp.n * sbp.n];
      
      status = mkl_sparse_d_mm(
        SPARSE_OPERATION_NON_TRANSPOSE, 1., F[i], da, 
        SPARSE_LAYOUT_COLUMN_MAJOR, MF[j],
        sbp.n, sbp.n * sbp.n, 
        1., ftmf, sbp.n);
      mkl_sparse_status(status);

      /*
      if (i * MF.size() + j == 8 || i * MF.size() + j == 21) {
        std::cout << i << " " << j << std::endl;
        std::cout << (i * MF.size() + j) * sbp.n * sbp.n <<  std::endl;
        for (std::size_t p = 0; p != sbp.n; ++p) {
          for (std::size_t k = 0; k != sbp.n; ++k) {
            std::cout << -ftmf[p * sbp.n + k]  << " ";
          }
          std::cout << std::endl;
        }
        std::cout << "mf " << j << std::endl;
        std::cout << std::endl;
        for (std::size_t p = 0; p != sbp.n * sbp.n; ++p) {
          for (std::size_t k = 0; k != sbp.n; ++k) {
            std::cout << MF[j][p * sbp.n + k]  << " ";
          }
          std::cout << std::endl;
        }
        
      }
      */
      /*
      for (std::size_t s = 0; s != sbp.n; ++s) {
        for (std::size_t k = 0; k != sbp.n; ++k) {
          std::cout << ftmf[s * sbp.n + k] << " "; // -= bc *-1 on D- FTMF
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      */
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

  /*
  int retval=PAPI_library_init(PAPI_VER_CURRENT);
  if (retval!=PAPI_VER_CURRENT) {
          fprintf(stderr,"Error initializing PAPI! %s\n",
                  PAPI_strerror(retval));
  }
  else {
      std::cout << "Loaded papi." << std::endl;
  }

  int eventset=PAPI_NULL;
  retval=PAPI_create_eventset(&eventset);
  if (retval!=PAPI_OK) {
    fprintf(stderr,"Error creating eventset! %s\n",
          PAPI_strerror(retval));
  }
  retval=PAPI_add_named_event(eventset,"PAPI_L3_TCM");
  if (retval!=PAPI_OK) {
    fprintf(stderr,"Error adding PAPI_L3_TCW: %s\n",
      PAPI_strerror(retval));
  }

  long long kount = 0;
  long long count = 0;
    
  

  PAPI_reset(eventset);
  retval=PAPI_start(eventset);
  if (retval!=PAPI_OK) {
    fprintf(stderr,"Error starting CUDA: %s\n",
      PAPI_strerror(retval));
  }
  */
  

  // Copy each block into a new array, sometimes wit
  std::size_t a, b, c, d, mindex, findex, ftindex, index;
  //#pragma omp parallel for private(a, b, c, d, mindex, findex, ftmf, lat) // num_threads(sbp.n_threads)
  for (std::size_t i = 0; i < interface_list.size(); i += 4) {
    a = interface_list[i];
    b = interface_list[i + 1];
    c = interface_list[i + 2];
    d = interface_list[i + 3];
    ftindex = FT_symbols[a][c] - 1;
    findex = F_symbols[c][b] - 1;
    auto r = c % sbp.n_blocks_dim == 0 ? 0 
          : c % sbp.n_blocks_dim == sbp.n_blocks_dim - 1 ? 2 
          : 1;
    // mindex = (r * 4) + F_symbols[c][b] - 1;
    // std::cout << i / 4 << " " << mindex << " " << findex << 
    //   " " << mindex * F.size() + findex <<std::endl;
    index = ((ftindex * MF.size()) + (r * 4) + findex) * sbp.n * sbp.n;
    ftmf = &FTMF[index];
    
    //std::cout << "FT " << FT_symbols[a][c] - 1 << std::endl;
    //std::cout << "F " << F_symbols[c][b] - 1 << std::endl;
    //std::cout << "M " << r << std::endl;
    //std::cout << a << " " << b << std::endl;

    /*
    if (a == 11 && b == 11) {
      for (std::size_t j = 0; j != sbp.n; ++j) {
        for (std::size_t k = 0; k != sbp.n; ++k) {
          std::cout << -ftmf[j * sbp.n + k]  << " ";
        }
        std::cout << std::endl;
      }
      
    }
    */
    //std::cout << index / (sbp.n * sbp.n) << std::endl;
    //std::cout << std::endl;

    for (std::size_t j = 0; j != sbp.n; ++j) {
      lat = &lambdaA[((j + (b * sbp.n)) * sbp.n_interfaces * sbp.n) + (a * sbp.n)];
      for (std::size_t k = 0; k != sbp.n; ++k) {
        lat[k] -= ftmf[j * sbp.n + k]; // -= bc *-1 on D- FTMF
        // kl = ((j + (a * sbp.n)) * sbp.n_interfaces * sbp.n) + (b * sbp.n) + k;
      }
    }
  }

  /*
  retval=PAPI_stop(eventset,&count);
  if (retval!=PAPI_OK) {
    fprintf(stderr,"Error stopping:  %s\n",
                          PAPI_strerror(retval));
  }
  else {
        printf("Measured %lld events \n",count);
        //printf("Measured %lld events \n",kount);
  }
  */
  

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
  // #pragma omp parallel for private(v) num_threads(sbp.n_threads)
  for (std::size_t j = 0; j != sbp.n_interfaces * sbp.n; ++j) {
    v = j % sbp.n;
    lambdaA[(j * sbp.n_interfaces * sbp.n) + j] += sbp.h1v[v] * 2. * sbp.τ;
    // lambdaA[(j * sbp.n_interfaces * sbp.n) + j] += 1;
  }
  

  mkl_free(FTMF);
  return;
}
